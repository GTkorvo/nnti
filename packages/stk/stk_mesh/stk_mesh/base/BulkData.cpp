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

#include <stk_mesh/base/BulkData.hpp>
#include <stddef.h>                     // for size_t, NULL
#include <string.h>                     // for memcpy, strcmp
#include <algorithm>                    // for sort, lower_bound, unique, etc
#include <boost/foreach.hpp>            // for auto_any_base, etc
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <sstream>
#include <fstream>
#include <iterator>                     // for back_insert_iterator, etc
#include <set>                          // for set, set<>::iterator, etc
#include <stk_mesh/base/Bucket.hpp>     // for Bucket, BucketIdComparator, etc
#include <stk_mesh/base/GetEntities.hpp>  // for get_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, print_entity_key, etc
#include <stk_mesh/baseImpl/EntityRepository.hpp>  // for EntityRepository, etc
#include <stk_mesh/baseImpl/Partition.hpp>  // for Partition
#include <stk_util/environment/ReportHandler.hpp>  // for ThrowRequireMsg, etc
#include <stk_util/parallel/ParallelComm.hpp>  // for CommBuffer, CommAll, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, all_reduce, etc
#include <stk_util/util/StaticAssert.hpp>  // for StaticAssert, etc
#include <stk_util/util/string_case_compare.hpp>
#include <string>                       // for char_traits, string, etc
#include <utility>                      // for pair, make_pair, swap
#include <vector>                       // for vector, etc
#include "boost/mpl/bool.hpp"           // for bool_
#include "boost/mpl/bool_fwd.hpp"       // for false_
#include "boost/unordered/detail/buckets.hpp"  // for iterator, etc
#include "stk_mesh/base/ConnectivityMap.hpp"  // for ConnectivityMap
#include "stk_mesh/base/Entity.hpp"     // for Entity, operator<<, etc
#include "stk_mesh/base/EntityCommDatabase.hpp"  // for pack_entity_info, etc
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, etc
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase, FieldMetaData, etc
#include "stk_mesh/base/FieldDataManager.hpp"  // for FieldDataManager, etc
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/Part.hpp"       // for Part, remove, etc
#include "stk_mesh/base/Relation.hpp"   // for Relation, etc
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Trace.hpp"      // for DiagIfWatching, Trace_, etc
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityRank, etc
#include "stk_mesh/baseImpl/BucketRepository.hpp"  // for BucketRepository
#include "stk_mesh/baseImpl/FieldRepository.hpp"  // for FieldVector
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/DistributedIndex.hpp"  // for DistributedIndex, etc
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine, etc
#include "stk_util/util/NamedPair.hpp"
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include "stk_util/util/SameType.hpp"   // for SameType, etc
#include "stk_util/util/TrackingAllocator.hpp"  // for tracking_allocator

namespace stk {
namespace mesh {

namespace impl {
int Counter::counter = 0;

#ifdef STK_MESH_ANALYZE_DYN_CONN
std::vector<DynConnData> DynConnMetrics::m_data;
#endif
}


///////////////////////////////////////////// Functions for creating entities

void fillEntityCommInfoForEntity(stk::mesh::Ghosting &ghost_id, stk::mesh::BulkData &mesh, std::vector<stk::mesh::EntityKey> nodes, EntityCommInfoVector &sharing_processors)
{
  size_t num_nodes = nodes.size();
  PairIterEntityComm initial_shared = mesh.entity_comm_map(nodes[0],ghost_id);
  sharing_processors.assign(initial_shared.first, initial_shared.second);

  for(size_t i = 1; i < num_nodes; ++i)
  {
    PairIterEntityComm tmp_shared  = mesh.entity_comm_map(nodes[i], ghost_id);
    EntityCommInfoVector new_shared_vec;

    std::set_intersection( sharing_processors.begin(), sharing_processors.end(),
                           tmp_shared.first, tmp_shared.second,
                           std::back_inserter(new_shared_vec) );
    sharing_processors = new_shared_vec;
    if (sharing_processors.empty())
    {
      break;
    }
  }
}

void fillSharedEntities(stk::mesh::Ghosting& ghost_id, stk::mesh::BulkData &mesh, std::vector<shared_entity_type> & shared_entity_map, std::vector<std::vector<shared_entity_type> > &shared_entities )
{
    for(std::vector<shared_entity_type>::const_iterator itr = shared_entity_map.begin(),
            end = shared_entity_map.end(); itr != end; ++itr)
    {
        EntityCommInfoVector sharing_processors;
        fillEntityCommInfoForEntity(ghost_id, mesh, itr->nodes, sharing_processors);

        for(EntityCommInfoVector::const_iterator comm_itr = sharing_processors.begin(),
                comm_end = sharing_processors.end(); comm_itr != comm_end; ++comm_itr)
        {
            if(comm_itr->proc != mesh.parallel_rank())
                shared_entities[comm_itr->proc].push_back(*itr);
        }
    }
}

void communicateSharedEntityInfo(stk::mesh::BulkData &mesh, CommAll &comm, std::vector<std::vector<shared_entity_type> > &shared_entities)
{
    for(int allocation_pass = 0; allocation_pass < 2; ++allocation_pass)
    {
        for(int proc = 0, parallel_size = mesh.parallel_size(); proc < parallel_size; ++proc)
        {
            if(proc != mesh.parallel_rank())
            {
                for(size_t e = 0, num_shared = shared_entities[proc].size(); e < num_shared; ++e)
                {
                    shared_entity_type const & sentity = shared_entities[proc][e];
                    size_t num_nodes_on_entity = sentity.nodes.size();
                    comm.send_buffer(proc).pack<stk::topology::topology_t>(sentity.topology);
                    //comm.send_buffer(proc).pack<size_t>(num_nodes_on_entity);  // TEMPORARY HACK: DO NOT PUSH
                    for (size_t i = 0; i < num_nodes_on_entity; ++i )
                    {
                        comm.send_buffer(proc).pack(sentity.nodes[i]);
                    }
                    comm.send_buffer(proc).pack(sentity.local_key);
                }
            }
        }

        if(allocation_pass == 0)
        {
            comm.allocate_buffers(mesh.parallel_size() / 4, 0);
        }
    }
    comm.communicate();
}

void unpackEntityInfromFromOtherProcsAndMarkEntitiesAsSharedAndTrackProcessorsThatNeedAlsoHaveEntity(CommAll &comm, stk::mesh::BulkData& mesh, std::vector<shared_entity_type> & shared_entity_map)
{
    for(int ip = mesh.parallel_size() - 1; ip >= 0; --ip)
    {
        if(ip != mesh.parallel_rank())
        {
            CommBuffer & buf = comm.recv_buffer(ip);
            while(buf.remaining())
            {
                shared_entity_type sentity;

                buf.unpack<stk::topology::topology_t>(sentity.topology);
                stk::topology entity_topology(sentity.topology);
                size_t num_nodes_on_entity = entity_topology.num_nodes();
                //size_t num_nodes_on_entity = 0;            // TEMPORARY HACK:
                //buf.unpack<size_t>(num_nodes_on_entity);   // DO NOT PUSH
                sentity.nodes.resize(num_nodes_on_entity);
                for (size_t i = 0; i < num_nodes_on_entity; ++i )
                {
                    buf.unpack<EntityKey>(sentity.nodes[i]);
                }
                buf.unpack<EntityKey>(sentity.global_key);

                std::vector<shared_entity_type>::iterator shared_itr = std::lower_bound(shared_entity_map.begin(), shared_entity_map.end(), sentity);

                //update the global global_key
                if(shared_itr != shared_entity_map.end() && *shared_itr == sentity)
                {
                    Entity entity = mesh.get_entity(shared_itr->local_key);
                    shared_itr->sharing_procs.push_back(ip);
                    if(ip < mesh.parallel_rank())
                    {
                        shared_itr->global_key = sentity.global_key;
                    }
                    mesh.mark_entity(entity, BulkData::IS_SHARED);
                }
            }
        }
    }
}

void BulkData::resolveUniqueIdForSharedEntityAndCreateCommMapInfoForSharingProcs(std::vector<shared_entity_type> & shared_entity_map)
{
   for(size_t i = 0, e = shared_entity_map.size(); i < e; ++i)
   {
       Entity entity = get_entity(shared_entity_map[i].local_key);
       if(shared_entity_map[i].global_key != shared_entity_map[i].local_key)
       {
           internal_change_entity_key(shared_entity_map[i].local_key, shared_entity_map[i].global_key, entity);
       }
       for(size_t j = 0; j < shared_entity_map[i].sharing_procs.size(); j++)
       {
           entity_comm_map_insert(entity, EntityCommInfo(stk::mesh::BulkData::SHARED, shared_entity_map[i].sharing_procs[j]));
       }
   }
}

void BulkData::update_shared_entities_global_ids(std::vector<shared_entity_type> & shared_entity_map)
{
    std::sort(shared_entity_map.begin(), shared_entity_map.end());

    // shared_edges[0] will contain all the edges this processor shares with processor 0
    std::vector<std::vector<shared_entity_type> > shared_entities(parallel_size());
    fillSharedEntities(shared_ghosting(), *this, shared_entity_map, shared_entities);

    CommAll comm(parallel());
    communicateSharedEntityInfo(*this, comm, shared_entities);
    unpackEntityInfromFromOtherProcsAndMarkEntitiesAsSharedAndTrackProcessorsThatNeedAlsoHaveEntity(comm, *this, shared_entity_map);
    resolveUniqueIdForSharedEntityAndCreateCommMapInfoForSharingProcs(shared_entity_map);
}
void BulkData::resolve_entity_sharing(stk::mesh::EntityRank entityRank, std::vector<EntityKey> &entity_keys)
{
    std::vector<shared_entity_type> shared_entities;
    impl::markEntitiesForResolvingSharingInfoUsingNodes(*this, entityRank, shared_entities);
    update_shared_entities_global_ids( shared_entities );

    for (size_t i=0; i<shared_entities.size();i++)
    {
        Entity entity = get_entity(shared_entities[i].global_key);
        if ( is_entity_marked(entity) == BulkData::IS_SHARED )
        {
            entity_keys.push_back(shared_entities[i].global_key);
        }
        mark_entity(entity, BulkData::NOT_MARKED);
    }
    std::sort(entity_keys.begin(), entity_keys.end());
}

/////////////////////////////////////// End functions for create edges

parallel::DistributedIndex::KeySpanVector
convert_entity_keys_to_spans( const MetaData & meta )
{
  // Make sure the distributed index can handle the EntityKey

  enum { OK = StaticAssert<
                SameType< uint64_t,
                          parallel::DistributedIndex::KeyType >::value >::OK };

  // Default constructed EntityKey has all bits set.

  const EntityKey invalid_key ;
  const EntityId  min_id = 1 ;
  const EntityId  max_id = invalid_key.id();

  const EntityRank rank_count = static_cast<EntityRank>(meta.entity_rank_count());

  parallel::DistributedIndex::KeySpanVector spans( rank_count );

  for ( EntityRank rank = stk::topology::NODE_RANK ; rank < rank_count ; ++rank ) {
    EntityKey key_min( rank , min_id );
    EntityKey key_max( rank , max_id );
    spans[rank].first  = key_min;
    spans[rank].second = key_max;
  }

  return spans ;
}

//----------------------------------------------------------------------

#ifdef STK_MESH_MODIFICATION_COUNTERS
unsigned BulkData::m_num_bulk_data_counter = 0;
#endif

BulkData::BulkData( MetaData & mesh_meta_data ,
                    ParallelMachine parallel
#ifdef SIERRA_MIGRATION
                    , bool add_fmwk_data
#endif
                    , ConnectivityMap const* arg_connectivity_map
                    , FieldDataManager *field_data_manager
                    , unsigned bucket_capacity
                    )
  :
#ifdef SIERRA_MIGRATION
    m_check_invalid_rels(true),
#endif
    m_entities_index( parallel, convert_entity_keys_to_spans(mesh_meta_data) ),
    m_entity_comm_map(),
    m_ghosting(),
    m_mesh_meta_data( mesh_meta_data ),
    m_parallel_machine( parallel ),
    m_parallel_size( parallel_machine_size( parallel ) ),
    m_parallel_rank( parallel_machine_rank( parallel ) ),
    m_sync_count( 0 ),
    m_sync_state( MODIFIABLE ),
    m_mark_entity(),
    m_add_node_sharing_called(false),
    m_closure_count(),
    m_entity_repo(*this),
    m_entity_comm_list(),
    m_volatile_fast_shared_comm_map(),
    m_ghost_parts(),
    m_deleted_entities(),
    m_deleted_entities_current_modification_cycle(),
    m_ghost_reuse_map(),
    m_meta_data_verified( false ),
    m_mesh_finalized(false),
    m_modification_begin_description("UNSET"),
    m_num_fields(-1), // meta data not necessarily committed yet
    m_keep_fields_updated(true),
    m_mesh_indexes(),
    m_entity_keys(),
    m_entity_states(),
    m_entity_sync_counts(),
    m_local_ids(),
#ifdef SIERRA_MIGRATION
    m_add_fmwk_data(add_fmwk_data),
    m_fmwk_bulk_ptr(NULL),
    m_fmwk_aux_relations(),
    m_fmwk_global_ids(),
    m_fmwk_shared_attrs(),
    m_fmwk_connect_counts(),
#endif
    m_default_field_data_manager(mesh_meta_data.entity_rank_count()),
    m_field_data_manager(field_data_manager),
    m_selector_to_buckets_map(),
    m_bucket_repository(
        *this,
        mesh_meta_data.entity_rank_count(),
        arg_connectivity_map != NULL ? *arg_connectivity_map :
        (mesh_meta_data.spatial_dimension() == 2 ? ConnectivityMap::default_map_2d() : ConnectivityMap::default_map()),
/*           (mesh_meta_data.spatial_dimension() == 2 ? ConnectivityMap::fixed_edges_map_2d() : ConnectivityMap::fixed_edges_map()) */
        bucket_capacity)
#ifdef STK_MESH_MODIFICATION_COUNTERS
    , m_num_bulk_data_counter++,
    m_modification_counters(),
    m_entity_modification_counters(),
#endif
#ifdef GATHER_GET_BUCKETS_METRICS
    , m_selector_to_count_map(),
    m_num_memoized_get_buckets_calls(0),
    m_num_non_memoized_get_buckets_calls(0),
    m_num_buckets_inserted_in_cache(0),
    m_num_buckets_removed_from_cache(0),
    m_num_modifications(0)
#endif
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
  std::ofstream outfile(create_modification_counts_filename().c_str());
  write_modification_labels_to_stream(outfile);
  outfile.close();
  reset_modification_counters();
#endif

  mesh_meta_data.set_mesh_bulk_data(this);

  if (m_field_data_manager == NULL)
  {
      m_field_data_manager = &m_default_field_data_manager;
  }

  initialize_arrays();

  m_ghost_parts.clear();
  internal_create_ghosting( "shared" );
  //shared part should reside in m_ghost_parts[0]
  internal_create_ghosting( "shared_aura" );

  m_sync_state = SYNCHRONIZED ;
}

void BulkData::reset_modification_counters()
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    for(unsigned j=0; j<static_cast<unsigned>(NumMethodTypes); ++j)
    {
        for(unsigned i=0; i<static_cast<unsigned>(NumModificationTypes); ++i)
        {
            m_modification_counters[j][i] = 0;
        }
        for(unsigned i=0; i<static_cast<unsigned>(NumEntityModificationTypes); ++i)
        {
            for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<mesh_meta_data().entity_rank_count(); rank++)
            {
                m_entity_modification_counters[j][rank][i] = 0;
            }
        }
    }
#endif
}

std::string BulkData::convert_label_for_method_type(const std::string &label, enum PublicOrInternalMethod methodType)
{
    std::string newLabel = label;
    if(methodType == INTERNAL)
    {
        newLabel = "INTERNAL-" + label;
    }
    return newLabel;
}

void BulkData::write_modification_entry_label(std::ostream& out, const std::string& label, enum PublicOrInternalMethod methodType)
{
    out << convert_label_for_method_type(label, methodType) << ", ";
}

void BulkData::write_entity_modification_entry_label(std::ostream& out, const std::string& label, enum PublicOrInternalMethod methodType)
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<mesh_meta_data().entity_rank_count(); rank++)
    {
        out << convert_label_for_method_type(label, methodType) <<"["<<rank<<"], ";
    }
    out << convert_label_for_method_type(label, methodType) << ", ";
#endif
}

void BulkData::write_modification_labels_to_stream_for_method_type(std::ostream& out, enum PublicOrInternalMethod methodType)
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    write_entity_modification_entry_label(out, "declare_entity", methodType);
    write_entity_modification_entry_label(out, "destroy_entity", methodType);
    write_entity_modification_entry_label(out, "change_entity_id", methodType);
    write_entity_modification_entry_label(out, "change_entity_parts", methodType);
    write_modification_entry_label(out, "change_entity_owner", methodType);
    write_modification_entry_label(out, "create_ghosting", methodType);
    write_modification_entry_label(out, "change_ghosting", methodType);
    write_modification_entry_label(out, "destroy_ghosting", methodType);
    write_modification_entry_label(out, "destroy_all_ghosting", methodType);
    write_modification_entry_label(out, "declare_relation", methodType);
    out << convert_label_for_method_type("destroy_relation", methodType);
#endif
}

void BulkData::write_modification_labels_to_stream(std::ostream& out)
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    out << "modification cycle, ";
    write_modification_labels_to_stream_for_method_type(out, PUBLIC);
    out << ", ";
    write_modification_labels_to_stream_for_method_type(out, INTERNAL);
    out << std::endl;
#endif
}

void BulkData::write_entity_modification_entry(std::ostream& out,
                                               enum PublicOrInternalMethod methodType,
                                               EntityModificationTypes entityModification)
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    unsigned sum = 0;
    for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<mesh_meta_data().entity_rank_count(); rank++)
    {
        out << m_entity_modification_counters[methodType][rank][entityModification]<<", ";
        sum += m_entity_modification_counters[methodType][rank][entityModification];
    }
    out << sum << ", ";
#endif
}

void BulkData::write_modification_counts_to_stream_for_method_type(std::ostream& out, enum PublicOrInternalMethod methodType)
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    write_entity_modification_entry(out, methodType, DECLARE_ENTITY);
    write_entity_modification_entry(out, methodType, DESTROY_ENTITY);
    write_entity_modification_entry(out, methodType, CHANGE_ENTITY_ID);
    write_entity_modification_entry(out, methodType, CHANGE_ENTITY_PARTS);
    out << m_modification_counters[methodType][CHANGE_ENTITY_OWNER]<<", ";
    out << m_modification_counters[methodType][CREATE_GHOSTING]<<", ";
    out << m_modification_counters[methodType][CHANGE_GHOSTING]<<", ";
    out << m_modification_counters[methodType][DESTROY_GHOSTING]<<", ";
    out << m_modification_counters[methodType][DESTROY_ALL_GHOSTING]<<", ";
    out << m_modification_counters[methodType][DECLARE_RELATION]<<", ";
    out << m_modification_counters[methodType][DESTROY_RELATION];
#endif
}

void BulkData::write_modification_counts_to_stream(std::ostream& out)
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    out << synchronized_count()<<", ";
    write_modification_counts_to_stream_for_method_type(out, PUBLIC);
    out << ", ";
    write_modification_counts_to_stream_for_method_type(out, INTERNAL);
    out << std::endl;
#endif
}

std::string BulkData::create_modification_counts_filename() const
{
    std::string fileName;
#ifdef STK_MESH_MODIFICATION_COUNTERS
    std::ostringstream oss;
    int numProcs = parallel_machine_size(MPI_COMM_WORLD);
    int procId = parallel_machine_rank(MPI_COMM_WORLD);
    oss<<"modification_counts_"<<m_num_bulk_data_counter<<"_np"<<numProcs<<"."<<procId<<".csv";
    fileName = oss.str();
#endif
    return fileName;
}

void BulkData::write_modification_counts()
{
#ifdef STK_MESH_MODIFICATION_COUNTERS
    std::ofstream outfile(create_modification_counts_filename().c_str(), std::ios::app);
    write_modification_counts_to_stream(outfile);
    reset_modification_counters();
#endif
}

#ifdef GATHER_GET_BUCKETS_METRICS
void BulkData::gather_and_print_get_buckets_metrics() const
{
  ParallelMachine world = MPI_COMM_WORLD; // HACK, but necessary to work with Fmwk
  CommAll comm_all(world);

  int real_rank = parallel_machine_rank(world);
  int real_size = parallel_machine_size(world);

  for (int phase = 0; phase < 2; ++phase) {
    CommBuffer& proc_0_buff = comm_all.send_buffer(0); // everything to rank 0

    proc_0_buff.pack<size_t>(m_num_buckets_inserted_in_cache);
    proc_0_buff.pack<size_t>(m_num_buckets_removed_from_cache);

    for (EntityRank r = 0; r < m_mesh_meta_data.entity_rank_count(); ++r) {
      proc_0_buff.pack<size_t>(buckets(r).size());
    }

    proc_0_buff.pack<size_t>(m_num_memoized_get_buckets_calls);
    proc_0_buff.pack<size_t>(m_num_non_memoized_get_buckets_calls);

    proc_0_buff.pack<size_t>(m_selector_to_count_map.size());
    for (SelectorCountMap::const_iterator itr = m_selector_to_count_map.begin(), end = m_selector_to_count_map.end(); itr != end; ++itr) {
      const EntityRank rank = itr->first.first;
      const size_t cache_hits = itr->second.first;
      const size_t bucket_trav_saved = itr->second.second;

      std::ostringstream out;
      out << itr->first.second;
      const std::string sel_str = out.str();

      proc_0_buff.pack<EntityRank>(rank);
      proc_0_buff.pack<size_t>(sel_str.size() + 1);
      proc_0_buff.pack<char>(sel_str.c_str(), sel_str.size() + 1);
      proc_0_buff.pack<size_t>(cache_hits);
      proc_0_buff.pack<size_t>(bucket_trav_saved);
    }

    if (phase == 0) { //allocation phase
      comm_all.allocate_buffers( m_parallel_size / 4 );
    }
    else { // communication phase
      comm_all.communicate();
    }
  }

  if (real_rank == 0) {
    std::ostringstream out;
    out << "get_buckets memoization_data:\n";
    out << "  Num modifications: " << m_num_modifications << "\n";

    typedef std::map< std::pair<EntityRank, std::string>, std::pair<size_t, size_t> > GlobalUsageMap;
    GlobalUsageMap global_usage_map;

    size_t global_num_memoized_get_buckets_calls_sum = 0;
    size_t global_num_non_memoized_get_buckets_calls_sum = 0;
    size_t global_num_cache_hits = 0;
    size_t global_num_bucket_trav_saved = 0;
    size_t global_num_buckets_inserted_in_cache = 0;
    size_t global_num_buckets_removed_from_cache = 0;

    std::vector<size_t> bucket_count(m_mesh_meta_data.entity_rank_count(), 0u);

    const size_t MAX_TEXT_LEN = 32768;
    char sel_text[ MAX_TEXT_LEN ];

    for ( int p = 0 ; p < real_size ; ++p ) {
      CommBuffer & buf = comm_all.recv_buffer(p);

      size_t num_buckets_inserted_in_cache = 0;
      size_t num_buckets_removed_from_cache = 0;

      buf.unpack<size_t>(num_buckets_inserted_in_cache);
      buf.unpack<size_t>(num_buckets_removed_from_cache);

      global_num_buckets_inserted_in_cache  += num_buckets_inserted_in_cache;
      global_num_buckets_removed_from_cache += num_buckets_removed_from_cache;

      for (EntityRank r = 0; r < m_mesh_meta_data.entity_rank_count(); ++r) {
        size_t count = 0;
        buf.unpack<size_t>(count);
        bucket_count[r] += count;
      }

      size_t num_memoized_get_buckets_calls = 0;
      size_t num_non_memoized_get_buckets_calls = 0;
      size_t map_size = 0;

      buf.unpack<size_t>(num_memoized_get_buckets_calls);
      buf.unpack<size_t>(num_non_memoized_get_buckets_calls);
      buf.unpack<size_t>(map_size);

      global_num_memoized_get_buckets_calls_sum     += num_memoized_get_buckets_calls;
      global_num_non_memoized_get_buckets_calls_sum += num_non_memoized_get_buckets_calls;

      for (size_t i = 0; i < map_size; ++i) {
        EntityRank rank = 0;
        size_t str_size = 0;
        size_t cache_hits = 0;
        size_t bucket_trav_saved = 0;

        buf.unpack<EntityRank>(rank);
        buf.unpack<size_t>(str_size);
        ThrowRequire(str_size < MAX_TEXT_LEN);
        buf.unpack<char>(sel_text, str_size);
        buf.unpack<size_t>(cache_hits);
        buf.unpack<size_t>(bucket_trav_saved);

        global_num_cache_hits += cache_hits;
        global_num_bucket_trav_saved += bucket_trav_saved;

        const std::string sel_str = sel_text;
        std::pair<EntityRank, std::string> search_key = std::make_pair(rank, sel_str);
        GlobalUsageMap::iterator f_itr = global_usage_map.find(search_key);
        if (f_itr == global_usage_map.end()) {
          global_usage_map[search_key] = std::make_pair(cache_hits, bucket_trav_saved);
        }
        else {
          f_itr->second.first  += cache_hits;
          f_itr->second.second += bucket_trav_saved;
        }
      }
    }

    size_t total_buckets = 0;
    out << "  Bucket counts by rank:\n";
    for (EntityRank r = 0; r < m_mesh_meta_data.entity_rank_count(); ++r) {
      out << "    " << r << ": " << bucket_count[r] << "\n";
      total_buckets += bucket_count[r];
    }
    out << "  Total buckets: " << total_buckets << "\n";
    out << "  Num buckets incrementally inserted in cache: " << global_num_buckets_inserted_in_cache << "\n";
    out << "  Num buckets incrementally removed from cache: " << global_num_buckets_removed_from_cache << "\n";
    out << "  Num memoized get_buckets calls: " << global_num_memoized_get_buckets_calls_sum << "\n";
    out << "  Num non-memoized get_buckets calls: " << global_num_non_memoized_get_buckets_calls_sum << "\n";
    out << "  Num cache hits: " << global_num_cache_hits << "\n";
    out << "  Num bucket traversals saved: " << global_num_bucket_trav_saved << "\n";
    out << "  Total number of selectors used: " << global_usage_map.size() << "\n";
    out << "  Selector usage counts:\n";
    for (GlobalUsageMap::const_iterator itr = global_usage_map.begin(), end = global_usage_map.end(); itr != end; ++itr) {
      out << "  (" << itr->first.first << ", " << itr->first.second << ") -> (" << itr->second.first << ", " << itr->second.second << ")\n";
    }
    std::cout << out.str() << std::endl;
  }
}
#endif

void BulkData::gather_and_print_mesh_partitioning() const
{
  ParallelMachine world = MPI_COMM_WORLD; // HACK, but necessary to work with Fmwk
  CommAll comm_all(world);

  int real_rank = parallel_machine_rank(world);
  int real_size = parallel_machine_size(world);

  for (int phase = 0; phase < 2; ++phase) {
    CommBuffer& proc_0_buff = comm_all.send_buffer(0); // everything to rank 0

    const EntityRank rank_count = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());

    for (EntityRank r = stk::topology::NODE_RANK; r < rank_count; ++r) {
      std::vector<impl::Partition*> partitions = m_bucket_repository.get_partitions(r);
      proc_0_buff.pack<size_t>(partitions.size());

      for (size_t p = 0, pe = partitions.size(); p < pe; ++p) {
        const std::vector<PartOrdinal> & parts = partitions[p]->get_legacy_partition_id();
        proc_0_buff.pack<size_t>(partitions[p]->size());

        std::ostringstream out;
        for (size_t i = 1, ie = parts.size() - 1; i < ie; ++i) {
          Part& part = m_mesh_meta_data.get_part(parts[i]);
          out << part.name() << "[" << part.primary_entity_rank() << "] ";
        }

        const std::string partition_str = out.str();

        proc_0_buff.pack<size_t>(partition_str.size() + 1);
        proc_0_buff.pack<char>(partition_str.c_str(), partition_str.size() + 1);
      }
    }

    if (phase == 0) { //allocation phase
      BABBLE_STK_PARALLEL_COMM(world, "          gather_and_print_mesh_partitioning calling allocate_buffers");
      comm_all.allocate_buffers( m_parallel_size / 4 );
    }
    else { // communication phase
      BABBLE_STK_PARALLEL_COMM(world, "          gather_and_pring_mesh_partitioning calling allocate_buffers");
      comm_all.communicate();
    }
  }

  if (real_rank == 0) {
    std::ostringstream out;
    out << "mesh partitioning data:\n";

    typedef std::map< std::pair<EntityRank, std::string>, size_t> GlobalPartitionMap;
    GlobalPartitionMap global_partition_map;

    const size_t MAX_TEXT_LEN = 32768;
    char partition_text[ MAX_TEXT_LEN ];

    for ( int p = 0 ; p < real_size ; ++p ) {
      CommBuffer & buf = comm_all.recv_buffer(p);

      const EntityRank rank_count = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());

      for (EntityRank r = stk::topology::NODE_RANK; r < rank_count; ++r) {
        size_t num_partitions_for_this_rank = 0;
        buf.unpack<size_t>(num_partitions_for_this_rank);

        for (size_t partition = 0; partition < num_partitions_for_this_rank; ++partition) {
          size_t num_entities_in_partition = 0;
          buf.unpack<size_t>(num_entities_in_partition);

          size_t str_size = 0;
          buf.unpack<size_t>(str_size);
          ThrowRequire(str_size < MAX_TEXT_LEN);
          buf.unpack<char>(partition_text, str_size);

          const std::string partition_str = partition_text;
          std::pair<EntityRank, std::string> search_key = std::make_pair(r, partition_str);
          GlobalPartitionMap::iterator f_itr = global_partition_map.find(search_key);
          if (f_itr == global_partition_map.end()) {
            global_partition_map[search_key] = num_entities_in_partition;
          }
          else {
            f_itr->second += num_entities_in_partition;
          }
        }
      }
    }

    for (GlobalPartitionMap::const_iterator itr = global_partition_map.begin(), end = global_partition_map.end(); itr != end; ++itr) {
      out << "  (" << itr->first.first << ", " << itr->first.second << ") -> " << itr->second << "\n";
    }
    std::cout << out.str() << std::endl;
  }
}

BulkData::~BulkData()
{

#ifdef STK_PROFILE_MEMORY
  ParallelMachine world = MPI_COMM_WORLD; // HACK, but necessary to work with Fmwk
  const int real_rank = parallel_machine_rank(world);
  print_max_stk_memory_usage(world, real_rank, std::cout);
#endif
#ifdef STK_MESH_ANALYZE_DYN_CONN
  print_dynamic_connectivity_profile(parallel(), parallel_rank(), std::cout);
#endif

#ifdef SIERRA_MIGRATION
  for(size_t i=0; i<m_fmwk_aux_relations.size(); ++i) {
    delete m_fmwk_aux_relations[i];
  }
#endif

#ifdef GATHER_GET_BUCKETS_METRICS
  gather_and_print_get_buckets_metrics();
#endif

#ifdef PRINT_MESH_PARTITIONING
  gather_and_print_mesh_partitioning();
#endif

  while ( ! m_ghosting.empty() ) {
    delete m_ghosting.back();
    m_ghosting.pop_back();
  }
}

void BulkData::get_selected_nodes(stk::mesh::Selector selector, stk::mesh::EntityVector& nodes)
{
    get_selected_entities( selector, this->buckets(stk::topology::NODE_RANK), nodes );
}

void BulkData::update_deleted_entities_container()
{
  //Question: should the m_deleted_entities container be sorted and uniqued?
  //I.e., should we guard against the same entity being deleted in consecutive modification cycles?

  while(!m_deleted_entities_current_modification_cycle.empty()) {
    size_t entity_offset = m_deleted_entities_current_modification_cycle.front();
    m_deleted_entities_current_modification_cycle.pop_front();
    m_deleted_entities.push_front(entity_offset);
  }

  // Reclaim offsets for deleted ghosted that were not regenerated
  for (GhostReuseMap::iterator m_itr = m_ghost_reuse_map.begin(), m_end = m_ghost_reuse_map.end(); m_itr != m_end; ++m_itr) {
    m_deleted_entities.push_front(m_itr->second);
  }

  m_ghost_reuse_map.clear();
}

size_t BulkData::total_field_data_footprint(EntityRank rank) const
{
  const std::vector< FieldBase * > & field_set = mesh_meta_data().get_fields();

  size_t retval = 0;
  for ( int i = 0; i < m_num_fields; ++i) {
    const FieldBase  & field = * field_set[i];
    retval += total_field_data_footprint(field, rank);
  }

  return retval;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::require_ok_to_modify() const
{
  ThrowRequireMsg( m_sync_state != SYNCHRONIZED,
                   "NOT in the ok-to-modify state" );
}

void BulkData::require_entity_owner( const Entity entity ,
                                     int owner ) const
{
  if (parallel_size() > 1 && bucket_ptr(entity) != NULL) {
    const bool error_not_owner = owner != parallel_owner_rank(entity) ;

    ThrowRequireMsg( !error_not_owner,
                     "Entity " << identifier(entity) << " owner is " <<
                     parallel_owner_rank(entity) << ", expected " << owner);
  }
}

void BulkData::require_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const
{
  const size_t rank_count = m_mesh_meta_data.entity_rank_count();
  const bool ok_id   = EntityKey::is_valid_id(ent_id);
  const bool ok_rank = ent_rank < rank_count && !(ent_rank == stk::topology::FACE_RANK && mesh_meta_data().spatial_dimension() == 2);

  ThrowRequireMsg( ok_rank,
                   "Bad key rank: " << ent_rank << " for id " << ent_id );

  ThrowRequireMsg( ok_id, "Bad id : " << ent_id);
}

bool BulkData::is_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const
{
  const size_t rank_count = m_mesh_meta_data.entity_rank_count();
  const bool ok_id   = EntityKey::is_valid_id(ent_id);
  const bool ok_rank = ent_rank < rank_count && !(ent_rank == stk::topology::FACE_RANK && mesh_meta_data().spatial_dimension() == 2);

  return ok_id && ok_rank;
}

void BulkData::require_metadata_committed()
{
  if (!m_mesh_meta_data.is_commit()) {
    m_mesh_meta_data.commit();
  }
}

//----------------------------------------------------------------------

const std::string& BulkData::get_modification_begin_description() const {
  return m_modification_begin_description;
}

//----------------------------------------------------------------------

bool BulkData::modification_begin(const std::string description)
{
  Trace_("stk::mesh::BulkData::modification_begin");

  parallel_machine_barrier( m_parallel_machine );

  ThrowRequireMsg( m_mesh_finalized == false, "Unable to modifiy, BulkData has been finalized.");

  if (m_sync_count == 0) {
    m_mesh_meta_data.set_mesh_on_fields(this);
  }

  if ( m_sync_state == MODIFIABLE && m_mesh_finalized == false ) return false ;


  m_modification_begin_description = description;


  if ( ! m_meta_data_verified ) {
    require_metadata_committed();

    if (parallel_size() > 1) {
      verify_parallel_consistency( m_mesh_meta_data , m_parallel_machine );
    }

    m_meta_data_verified = true ;
  }
  else {
    ++m_sync_count ;

    for(unsigned i=0, iend=m_entity_states.size(); i<iend; ++i) {
      if(m_entity_states[i] != Deleted) {
        m_entity_states[i] = Unchanged;
      }
    }
  }

  // // It might be overkill to call this on every modification cycle.

  m_sync_state = MODIFIABLE ;

  return true ;
}

namespace impl {

struct MarkAsModified
{
    MarkAsModified(BulkData & mesh_in) : mesh(mesh_in) {}
    void operator()(Entity entity){
        mesh.set_state(entity, Modified);
    }
    BulkData & mesh;
};

struct OnlyVisitUnchanged
{
    OnlyVisitUnchanged(BulkData & mesh_in) : mesh(mesh_in) {}
    bool operator()(Entity entity){
        if (mesh.state(entity) == Unchanged) {
            return true;
        }
        return false;
    }
    BulkData & mesh;
};

} //namespace impl

void BulkData::mark_entity_and_upward_related_entities_as_modified(Entity entity)
{
  TraceIfWatching("stk::mesh::BulkData::log_modified_and_propagate", LOG_ENTITY, entity_key(entity));

  impl::MarkAsModified mam(*this);
  impl::OnlyVisitUnchanged ovu(*this);
  impl::VisitUpwardClosureGeneral(*this, entity, mam, ovu);
}

size_t BulkData::count_relations(Entity entity) const
{
  const MeshIndex &mesh_idx = mesh_index(entity);

  const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());
  size_t count = 0;
  for (EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
  {
    count += mesh_idx.bucket->num_connectivity(mesh_idx.bucket_ordinal, irank);
  }
  return count;
}

bool BulkData::has_no_relations(Entity entity) const
{
  const MeshIndex &mesh_idx = mesh_index(entity);

  const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());
  for (EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
  {
    if (mesh_idx.bucket->num_connectivity(mesh_idx.bucket_ordinal, irank) > 0)
    {
      return false;
    }
  }
  return true;
}

unsigned BulkData::count_valid_connectivity(Entity entity, EntityRank rank) const
{
  if (bucket(entity).connectivity_type(rank) == FIXED_CONNECTIVITY) {

    m_check_invalid_rels = false;
    Entity const *rel_iter = begin(entity, rank);
    Entity const *rel_end = end(entity, rank);
    m_check_invalid_rels = true;

    unsigned count = 0;
    for (; rel_iter != rel_end; ++rel_iter)
    {
      if (rel_iter->is_local_offset_valid())
      {
        ++count;
      }
    }
    return count;
  }
  else {
    return bucket(entity).num_connectivity(bucket_ordinal(entity), rank);
  }
}

unsigned BulkData::count_valid_connectivity(Entity entity) const
{
  unsigned count = 0;
  const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());
  for (EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
  {
    count += count_valid_connectivity(entity, irank);
  }
  return count;
}

size_t BulkData::generate_next_local_offset(size_t preferred_offset)
{
  size_t new_local_offset = m_mesh_indexes.size();

  if (preferred_offset != 0) {
    new_local_offset = preferred_offset;
  }
  else if (!m_deleted_entities.empty()) {
    new_local_offset = m_deleted_entities.front();
    m_deleted_entities.pop_front();
  }

  MeshIndex mesh_index = {NULL, 0};
  EntityKey invalid_key;

  if (new_local_offset == m_mesh_indexes.size()) {
    m_mesh_indexes.push_back(mesh_index);
    m_entity_keys.push_back(invalid_key);
    m_entity_states.push_back(Created);
    m_mark_entity.push_back(NOT_MARKED);
    m_closure_count.push_back(static_cast<uint16_t>(0));
    m_entity_sync_counts.push_back(0);
    m_local_ids.push_back(stk::mesh::GetInvalidLocalId());

#ifdef SIERRA_MIGRATION
    if (m_add_fmwk_data) {
      m_fmwk_aux_relations.push_back(NULL);
      m_fmwk_global_ids.push_back(0);
      m_fmwk_shared_attrs.push_back(NULL);
      m_fmwk_connect_counts.push_back(0);
    }
#endif
  }
  else {
    //re-claiming space from a previously-deleted entity:

    m_mesh_indexes[new_local_offset] = mesh_index;
    m_entity_keys[new_local_offset] = invalid_key;
    m_mark_entity[new_local_offset] = NOT_MARKED;
    m_entity_states[new_local_offset] = Created;
    m_closure_count[new_local_offset] = static_cast<uint16_t>(0);
    m_entity_sync_counts[new_local_offset] = 0;
    m_local_ids[new_local_offset] = stk::mesh::GetInvalidLocalId();

#ifdef SIERRA_MIGRATION
    if (m_add_fmwk_data) {
      //bulk-data allocated aux-relation vector, so delete it here.
      delete m_fmwk_aux_relations[new_local_offset];
      m_fmwk_aux_relations[new_local_offset] = NULL;
      m_fmwk_global_ids[new_local_offset] = 0;
      //don't delete shared-attr, it was allocated by fmwk.
      m_fmwk_shared_attrs[new_local_offset] = NULL;
      m_fmwk_connect_counts[new_local_offset] = 0;
    }
#endif
  }

  return new_local_offset;
}

void BulkData::initialize_arrays()
{
  ThrowRequireMsg((m_mesh_indexes.size() == 0) && (m_entity_keys.size() == 0)
                        && (m_entity_states.size() == 0) && (m_entity_sync_counts.size() == 0),
                   "BulkData::initialize_arrays() called by something other than constructor");

  MeshIndex mesh_index = {NULL, 0};
  m_mesh_indexes.push_back(mesh_index);

  EntityKey invalid_key;
  m_entity_keys.push_back(invalid_key);

  m_entity_states.push_back(Deleted);
  m_mark_entity.push_back(NOT_MARKED);
  m_closure_count.push_back(static_cast<uint16_t>(0));
  m_entity_sync_counts.push_back(0);
  m_local_ids.push_back(stk::mesh::GetInvalidLocalId());

#ifdef SIERRA_MIGRATION
  if (m_add_fmwk_data) {
    m_fmwk_aux_relations.push_back(NULL);
    m_fmwk_global_ids.push_back(0);
    m_fmwk_shared_attrs.push_back(NULL);
    m_fmwk_connect_counts.push_back(0);
  }
#endif
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// The add_parts must be full ordered and consistent,
// i.e. no bad parts, all supersets included, and
// owner & used parts match the owner value.

//----------------------------------------------------------------------

Entity BulkData::declare_entity( EntityRank ent_rank , EntityId ent_id)
{
    INCREMENT_ENTITY_MODIFICATION_COUNTER(PUBLIC, ent_rank, DECLARE_ENTITY);
    PartVector parts(1, &mesh_meta_data().universal_part());
    return internal_declare_entity(ent_rank, ent_id, parts);
}

Entity BulkData::declare_entity( EntityRank ent_rank , EntityId ent_id , Part& part)
{
    INCREMENT_ENTITY_MODIFICATION_COUNTER(PUBLIC, ent_rank, DECLARE_ENTITY);
    PartVector parts(1, &part);
    return internal_declare_entity( ent_rank, ent_id, parts);
}

Entity BulkData::declare_entity( EntityRank ent_rank , EntityId ent_id ,
                                 const PartVector & parts )
{
    INCREMENT_ENTITY_MODIFICATION_COUNTER(PUBLIC, ent_rank, DECLARE_ENTITY);
    return internal_declare_entity(ent_rank, ent_id, parts);
}

void BulkData::add_node_sharing( Entity node, int sharing_proc )
{
  // Only valid to specify sharing information for newly-created nodes
  ThrowRequire(entity_rank(node) == stk::topology::NODE_RANK);
  ThrowRequire(state(node) == Created);

  m_add_node_sharing_called = true;

  mark_entity(node, IS_SHARED);
  entity_comm_map_insert(node, EntityCommInfo(stk::mesh::BulkData::SHARED, sharing_proc));
}

void BulkData::internal_add_node_sharing( Entity node, int sharing_proc )
{
  // Only valid to specify sharing information for newly-created nodes
  ThrowRequire(entity_rank(node) == stk::topology::NODE_RANK);

  m_add_node_sharing_called = true;

  mark_entity(node, IS_SHARED);
  entity_comm_map_insert(node, EntityCommInfo(stk::mesh::BulkData::SHARED, sharing_proc));
}

void BulkData::internal_remove_node_sharing( EntityKey node, int sharing_proc )
{
    bool result = entity_comm_map_erase(node, EntityCommInfo(stk::mesh::BulkData::SHARED, sharing_proc));
    ThrowRequireMsg(result, "ERROR, entity_comm_map_erase called for node " << node << " and proc " << sharing_proc << " and returned false" );
//    PairIterEntityComm piec = entity_comm_map_shared(node);
//    if (piec.begin() == piec.end()) {
//        mark_entity(get_entity(node), NOT_MARKED);
//    }
}

Entity BulkData::internal_declare_entity( EntityRank ent_rank , EntityId ent_id ,
                                 const PartVector & parts )
{
    m_check_invalid_rels = false;

  require_ok_to_modify();

  require_good_rank_and_id(ent_rank, ent_id);

  EntityKey key( ent_rank , ent_id );
  TraceIfWatching("stk::mesh::BulkData::internal_declare_entity", LOG_ENTITY, key);
  DiagIfWatching(LOG_ENTITY, key, "declaring entity with parts " << parts);

  std::pair< Entity , bool > result = internal_create_entity( key );

  Entity declared_entity = result.first;

  if ( !result.second) {
    // An existing entity, the owner must match.
    require_entity_owner( declared_entity , m_parallel_rank );
    DiagIfWatching(LOG_ENTITY, key, "existing entity: " << entity_key(declared_entity));
  }

  //------------------------------

  Part * const owns = & m_mesh_meta_data.locally_owned_part();

  PartVector rem ;
  PartVector add( parts );
  add.push_back( owns );

  internal_verify_and_change_entity_parts( declared_entity , add , rem );

  if ( result.second ) {
    this->set_parallel_owner_rank_but_not_comm_lists(declared_entity, m_parallel_rank);
    set_synchronized_count(declared_entity, m_sync_count);
    DiagIfWatching(LOG_ENTITY, key, "new entity: " << entity_key(declared_entity));
  }

  m_check_invalid_rels = true;

  return declared_entity ;
}

void BulkData::change_entity_id( EntityId id, Entity entity)
{
// When stk parallel is used within Fmwk, this assertion is violated
#ifndef SIERRA_MIGRATION
  ThrowAssertMsg(parallel_size() == 1,
                 "change_entity_id only supported in serial");
#endif

  EntityRank e_rank = entity_rank(entity);
  INCREMENT_ENTITY_MODIFICATION_COUNTER(PUBLIC, e_rank, CHANGE_ENTITY_ID);

  require_ok_to_modify();
  require_good_rank_and_id(e_rank, id);

  EntityKey new_key(e_rank,id);
  EntityKey old_key = entity_key(entity);

  internal_change_entity_key(old_key, new_key, entity);
}

void BulkData::internal_change_entity_key( EntityKey old_key, EntityKey new_key, Entity entity)
{
  INCREMENT_ENTITY_MODIFICATION_COUNTER(INTERNAL, old_key.rank(), CHANGE_ENTITY_ID);

  m_entity_repo.update_entity_key(new_key, old_key, entity);
  set_entity_key(entity, new_key);
  this->bucket(entity).getPartition()->set_flag_needs_to_be_sorted(true);
}

//----------------------------------------------------------------------

bool BulkData::destroy_entity( Entity entity, bool was_ghost )
{
    INCREMENT_ENTITY_MODIFICATION_COUNTER(PUBLIC, entity_rank(entity), DESTROY_ENTITY);
    return internal_destroy_entity(entity, was_ghost);
}

bool BulkData::internal_destroy_entity( Entity entity, bool was_ghost )
{
  const EntityKey key = entity_key(entity);
  TraceIfWatching("stk::mesh::BulkData::destroy_entity", LOG_ENTITY, key);
  DiagIfWatching(LOG_ENTITY, key, "entity state: " << key);

  require_ok_to_modify();

  m_check_invalid_rels = false;

  if (!is_valid(entity)) {
    m_check_invalid_rels = true;
    return false;
  }

  const bool ghost = was_ghost || in_receive_ghost(key);
  const EntityRank erank = entity_rank(entity);

  INCREMENT_ENTITY_MODIFICATION_COUNTER(INTERNAL, erank, DESTROY_ENTITY);

  const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());
  for (EntityRank irank = static_cast<EntityRank>(erank + 1); irank != end_rank; ++irank) {
    if (num_connectivity(entity, irank) > 0) {
      m_check_invalid_rels = true;
      return false;
    }
  }

  //------------------------------
  // Immediately remove it from relations and buckets.
  // Postpone deletion until modification_end to be sure that
  // 1) No attempt is made to re-create it.
  // 2) Parallel index is cleaned up.
  // 3) Parallel sharing is cleaned up.
  // 4) Parallel ghosting is cleaned up.
  //
  // Must clean up the parallel lists before fully deleting the entity.

  // It is important that relations be destroyed from highest to lowest rank so
  // that the back relations are destroyed first.
  EntityVector temp_entities;
  std::vector<ConnectivityOrdinal> temp_ordinals;
  Entity const* rel_entities = NULL;
  int num_conn = 0;
  ConnectivityOrdinal const* rel_ordinals;
  for (EntityRank irank = end_rank; irank != stk::topology::BEGIN_RANK; )
  {
    --irank;

    if (connectivity_map().valid(erank, irank)) {
      num_conn     = num_connectivity(entity, irank);
      rel_entities = begin(entity, irank);
      rel_ordinals = begin_ordinals(entity, irank);
    }
    else {
      num_conn     = get_connectivity(*this, entity, irank, temp_entities, temp_ordinals);
      rel_entities = &*temp_entities.begin();
      rel_ordinals = &*temp_ordinals.begin();
    }

    for (int j = num_conn; j > 0; )
    {
      --j;
      if (is_valid(rel_entities[j])) {
        internal_destroy_relation(entity, rel_entities[j], rel_ordinals[j]);
      }
    }
  }

  // If this is a ghosted entity, store key->local_offset so that local_offset can be
  // reused if the entity is recreated in the next aura-regen. This will prevent clients
  // from having their handles to ghosted entities go invalid when the ghost is refreshed.
  if ( ghost ) {
    DiagIfWatching(LOG_ENTITY, key, "Adding to ghost reuse map");
    m_ghost_reuse_map[key] = entity.local_offset();
  }

  // We need to save these items and call remove_entity AFTER the call to
  // destroy_later because remove_entity may destroy the bucket
  // which would cause problems in m_entity_repo.destroy_later because it
  // makes references to the entity's original bucket.

  // Need to invalidate Entity handles in comm-list
  EntityCommListInfoVector::iterator lb_itr =
    std::lower_bound(m_entity_comm_list.begin(), m_entity_comm_list.end(), key);
  if (lb_itr != m_entity_comm_list.end() && lb_itr->key == key) {
    lb_itr->entity = Entity();
  }

  remove_entity_callback(erank, bucket(entity).bucket_id(), bucket_ordinal(entity));

  m_entities_index.register_removed_key( key );

  bucket(entity).getPartition()->remove(entity);
  m_entity_repo.destroy_entity(key, entity );
  m_entity_states[entity.local_offset()] = Deleted;
  m_mark_entity[entity.local_offset()] = NOT_MARKED;
  m_closure_count[entity.local_offset()] = static_cast<uint16_t>(0u);
  if ( !ghost ) {
    m_deleted_entities_current_modification_cycle.push_front(entity.local_offset());
  }

  m_check_invalid_rels = true;
  return true ;
}

//----------------------------------------------------------------------

void BulkData::generate_new_entities(const std::vector<size_t>& requests,
                                 std::vector<Entity>& requested_entities)
// requests = number of nodes needed, number of elements needed, etc.
{
    Trace_("stk::mesh::BulkData::generate_new_entities");

    typedef stk::parallel::DistributedIndex::KeyType KeyType;
    typedef stk::parallel::DistributedIndex::KeyTypeVector KeyTypeVector;
    typedef std::vector<KeyTypeVector> RequestKeyVector;

    RequestKeyVector requested_key_types;

    m_entities_index.generate_new_keys(requests, requested_key_types);

    //generating 'owned' entities
    Part * const owns = &m_mesh_meta_data.locally_owned_part();

    std::vector<Part*> rem;
    std::vector<Part*> add;
    add.push_back(owns);

    requested_entities.clear();
    unsigned total_number_of_key_types_requested = 0;
    for(RequestKeyVector::const_iterator itr = requested_key_types.begin(); itr != requested_key_types.end(); ++itr)
    {
        const KeyTypeVector& key_types = *itr;
        for(KeyTypeVector::const_iterator kitr = key_types.begin(); kitr != key_types.end(); ++kitr)
        {
            ++total_number_of_key_types_requested;
        }
    }
    requested_entities.reserve(total_number_of_key_types_requested);
    addMeshEntities(requested_key_types, rem, add, requested_entities);
}

std::pair<Entity, bool> BulkData::internal_create_entity(EntityKey key, size_t preferred_offset)
{
    INCREMENT_ENTITY_MODIFICATION_COUNTER(INTERNAL, key.rank(), DECLARE_ENTITY);
    return m_entity_repo.internal_create_entity(key, preferred_offset);
}

void BulkData::addMeshEntities(const std::vector< stk::parallel::DistributedIndex::KeyTypeVector >& requested_key_types,
       const std::vector<Part*> &rem, const std::vector<Part*> &add, std::vector<Entity>& requested_entities)
{
    typedef stk::parallel::DistributedIndex::KeyType       KeyType;
    typedef stk::parallel::DistributedIndex::KeyTypeVector KeyTypeVector;
    typedef std::vector< KeyTypeVector > RequestKeyVector;

    for(RequestKeyVector::const_iterator itr = requested_key_types.begin(); itr != requested_key_types.end(); ++itr)
    {
        const KeyTypeVector & key_types = *itr;
        for(KeyTypeVector::const_iterator kitr = key_types.begin(); kitr != key_types.end(); ++kitr)
        {
            EntityKey key(static_cast<EntityKey::entity_key_t>((*kitr)));
            require_good_rank_and_id(key.rank(), key.id());
            std::pair<Entity, bool> result = internal_create_entity(key);

            //if an entity is declared with the declare_entity function in
            //the same modification cycle as the generate_new_entities
            //function, and it happens to generate a key that was declared
            //previously in the same cycle it is an error
            ThrowErrorMsgIf( ! result.second,
                    "Generated id " << key.id() << " of rank " << key.rank() <<
                    " which was already used in this modification cycle.");

            // A new application-created entity

            Entity new_entity = result.first;

            //add entity to 'owned' part
            internal_verify_and_change_entity_parts(new_entity, add, rem);
            requested_entities.push_back(new_entity);

            this->set_parallel_owner_rank_but_not_comm_lists(new_entity, m_parallel_rank);
            set_synchronized_count(new_entity, m_sync_count);
        }
    }
}


bool BulkData::in_shared(EntityKey key, int proc) const
{
  PairIterEntityComm sharing = entity_comm_map_shared(key);
  for ( ; !sharing.empty(); ++sharing ) {
    if ( proc == sharing->proc ) {
      return true ;
    }
  }
  return false ;
}

bool BulkData::is_ghosted_onto_proc( EntityKey key, int otherProc ) const
{
  const int proc = parallel_rank();
  const int owner_rank = entity_comm_map_owner(key);
  if ( proc == owner_rank )
  {
      for ( PairIterEntityComm ec = entity_comm_map(key); ! ec.empty() ; ++ec ) {
        if ( ec->ghost_id == 1 &&
             ec->proc  == otherProc &&
             ec->proc     != proc ) {
          return true;
        }
      }
  }
  return false;
}

bool BulkData::is_ghosted_onto_another_proc( EntityKey key ) const
{
  const int proc = parallel_rank();
  const int owner_rank = entity_comm_map_owner(key);
  if ( proc == owner_rank )
  {
      for ( PairIterEntityComm ec = entity_comm_map(key); ! ec.empty() ; ++ec ) {
        if ( ec->ghost_id == 1 &&
             ec->proc     != proc ) {
          return true;
        }
      }
  }
  return false;
}

bool BulkData::in_send_ghost( EntityKey key , int proc ) const
{
  const int owner_rank = entity_comm_map_owner(key);
  for ( PairIterEntityComm ec = entity_comm_map(key); ! ec.empty() ; ++ec ) {
    if ( ec->ghost_id != 0 &&
         ec->proc     != owner_rank &&
         ec->proc     == proc ) {
      return true;
    }
  }
  return false;
}

bool BulkData::in_ghost( const Ghosting & ghost , EntityKey key , int proc ) const
{
  // Ghost communication from owner.
  EntityCommInfo tmp( ghost.ordinal() , proc );

  PairIterEntityComm ec = entity_comm_map(key);
  EntityCommInfoVector::const_iterator i =
    std::lower_bound( ec.begin(), ec.end() , tmp );

  return i != ec.end() && tmp == *i ;
}

void BulkData::comm_procs( EntityKey key, std::vector<int> & procs ) const
{
  procs.clear();
  for ( PairIterEntityComm ec = entity_comm_map(key); ! ec.empty() ; ++ec ) {
    procs.push_back( ec->proc );
  }
  std::sort( procs.begin() , procs.end() );
  std::vector<int>::iterator
    i = std::unique( procs.begin() , procs.end() );
  procs.erase( i , procs.end() );
}
void BulkData::comm_shared_procs( EntityKey key, std::vector<int> & procs ) const
{
  procs.clear();
  for ( PairIterEntityComm ec = entity_comm_map_shared(key); ! ec.empty() ; ++ec ) {
    procs.push_back( ec->proc );
  }
  std::sort( procs.begin() , procs.end() );
  std::vector<int>::iterator
    i = std::unique( procs.begin() , procs.end() );
  procs.erase( i , procs.end() );
}

void BulkData::shared_procs_intersection( std::vector<EntityKey> & keys, std::vector<int> & procs ) const
{

  procs.clear();
  int num = keys.size();
  std::vector<int> procs_tmp;
  for (int i = 0; i < num; ++i)
  {

    comm_shared_procs(keys[i], procs_tmp);

    if (i == 0)
      procs.swap(procs_tmp);
    else
    {
      // subsequent loops keep the intersection
      std::vector<int> result;
      std::back_insert_iterator<std::vector<int> > result_itr(result);
      std::set_intersection(procs.begin(),
                            procs.end(),
                            procs_tmp.begin(),
                            procs_tmp.end(),
                            result_itr,
                            std::less<int>());
      procs.swap(result);
    }
  }
}

void BulkData::comm_procs( const Ghosting & ghost ,
                           EntityKey key, std::vector<int> & procs ) const
{
  procs.clear();
  for ( PairIterEntityComm ec = entity_comm_map(key); ! ec.empty() ; ++ec ) {
    if ( ec->ghost_id == ghost.ordinal() ) {
      procs.push_back( ec->proc );
    }
  }
}

void BulkData::internal_change_owner_in_comm_data(const EntityKey& key, int new_owner)
{
  const bool changed = m_entity_comm_map.change_owner_rank(key, new_owner);
  if (changed) {
    EntityCommListInfoVector::iterator lb_itr = std::lower_bound(m_entity_comm_list.begin(),
                                                                        m_entity_comm_list.end(),
                                                                        key);
    if (lb_itr != m_entity_comm_list.end() && lb_itr->key == key) {
      lb_itr->owner = new_owner;
    }
  }
}

void BulkData::internal_sync_comm_list_owners()
{
  for (size_t i = 0, e = m_entity_comm_list.size(); i < e; ++i) {
    m_entity_comm_list[i].owner = parallel_owner_rank(m_entity_comm_list[i].entity);
  }
}

void BulkData::deactivate_field_updating()
{
  if (m_num_fields > -1) {
    //if fields have already been allocated, then we can't deactivate the updating
    //of field-data.
    m_keep_fields_updated = true;
    return;
  }

  m_keep_fields_updated = false;
}

void BulkData::allocate_field_data()
{
  if (m_keep_fields_updated == true) {
    //fields are already allocated, nothing to do here.
    return;
  }

  //temporary (hopefully) kludge:
  //calling the buckets(rank) getter causes partitions/buckets to potentially
  //be reorganized (including deleting buckets) and so we need to do it
  //before flipping the m_keep_fields_updated flag...
  for(EntityRank rank = stk::topology::NODE_RANK; rank < mesh_meta_data().entity_rank_count(); ++rank) {
    this->buckets(rank);
  }

  m_keep_fields_updated = true;
  //now loop over all buckets and call the 'new_bucket_callback' method which
  //will allocate field-data for that bucket.

  const std::vector< FieldBase * > & field_set = mesh_meta_data().get_fields();
  if (m_num_fields == -1) {
    // hasn't been set yet
    m_num_fields = field_set.size();
  }

  for(EntityRank rank = stk::topology::NODE_RANK; rank < mesh_meta_data().entity_rank_count(); ++rank) {
      const std::vector<Bucket*>& buckets = this->buckets(rank);
      m_field_data_manager->allocate_field_data(rank, buckets, field_set);
  }
}

void BulkData::new_bucket_caching(EntityRank rank, Bucket* new_bucket)
{
    // update selector map
    if (new_bucket != NULL) {
      for (SelectorBucketMap::iterator itr = m_selector_to_buckets_map.begin(), end = m_selector_to_buckets_map.end();
           itr != end; ++itr) {
        Selector const& sel = itr->first.second;
        const EntityRank map_rank = itr->first.first;
        if (map_rank == rank && sel(*new_bucket)) {
  #ifdef GATHER_GET_BUCKETS_METRICS
          ++m_num_buckets_inserted_in_cache;
  #endif
          TrackedBucketVector& cached_buckets = itr->second;
          TrackedBucketVector::iterator lb_itr = std::lower_bound(cached_buckets.begin(), cached_buckets.end(), new_bucket, BucketIdComparator());
          cached_buckets.insert(lb_itr, new_bucket);
        }
      }
    }
}

void BulkData::new_bucket_callback(EntityRank rank, const PartVector& superset_parts, size_t capacity, Bucket* new_bucket)
{
  this->new_bucket_caching(rank, new_bucket);

  if (!m_keep_fields_updated) {
    return;
  }

  const std::vector< FieldBase * > & field_set = mesh_meta_data().get_fields();

  if (m_num_fields == -1) {
    // hasn't been set yet
    m_num_fields = field_set.size();
  }

  m_field_data_manager->allocate_bucket_field_data(rank, field_set, superset_parts, capacity);
}

//
//  Copy fields from src to dst entity.  If the field size of either entity is zero, do nothing.  If the field
//  size of of both entities are non-zero, then the sizes must match
//

void BulkData::copy_entity_fields_callback(EntityRank dst_rank, unsigned dst_bucket_id, Bucket::size_type dst_bucket_ord,
                                           unsigned src_bucket_id, Bucket::size_type src_bucket_ord, const std::vector<FieldBase*>* field_set)
{
    //
    //  If field set is passed in copy only the defined fields.  Also assume the fields are valid for the bucket
    //
    if(field_set)
    {
        for(int i = 0, iend = field_set->size(); i < iend; ++i)
        {
            const int src_size = (*field_set)[i]->get_meta_data_for_field()[src_bucket_id].m_bytes_per_entity;
            unsigned char * const src = (*field_set)[i]->get_meta_data_for_field()[src_bucket_id].m_data;
            unsigned char * const dst = (*field_set)[i]->get_meta_data_for_field()[dst_bucket_id].m_data;

            ThrowAssert(src_size == (*field_set)[i]->get_meta_data_for_field()[dst_bucket_id].m_bytes_per_entity);

            std::memcpy(dst + src_size * dst_bucket_ord,
                    src + src_size * src_bucket_ord,
                    src_size);
        }
    }
    else
    {
        if(!m_keep_fields_updated)
        {
            return;
        }

        const std::vector<FieldBase *>& allFields = mesh_meta_data().get_fields((stk::topology::rank_t) dst_rank);
        for(int i = 0, iend = allFields.size(); i < iend; ++i)
        {
            const int src_size = allFields[i]->get_meta_data_for_field()[src_bucket_id].m_bytes_per_entity;
            if(src_size == 0)
            {
                continue;
            }

            unsigned char * const src = allFields[i]->get_meta_data_for_field()[src_bucket_id].m_data;
            const int dst_size = allFields[i]->get_meta_data_for_field()[dst_bucket_id].m_bytes_per_entity;

            if(dst_size)
            {
                unsigned char * const dst = allFields[i]->get_meta_data_for_field()[dst_bucket_id].m_data;
                ThrowAssertMsg( dst_size == src_size,
                        "Incompatible field sizes: " << dst_size << " != " << src_size);

                std::memcpy(dst + dst_size * dst_bucket_ord,
                        src + src_size * src_bucket_ord,
                        dst_size);
            }
        }
    }
}

void BulkData::add_entity_callback(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord)
{
    if(!m_keep_fields_updated)
    {
        return;
    }
    const std::vector<FieldBase *> &fields = mesh_meta_data().get_fields();
    m_field_data_manager->add_field_data_for_entity(fields, rank, bucket_id, bucket_ord);
}

void BulkData::remove_entity_field_data_callback(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord)
{
    if (!m_keep_fields_updated) {
      return;
    }
    const std::vector<FieldBase *> &fields = mesh_meta_data().get_fields();
    m_field_data_manager->remove_field_data_for_entity(rank, bucket_id, bucket_ord, fields);
}

void BulkData::remove_entity_callback(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord)
{
  if (!m_keep_fields_updated) {
    return;
  }
  const std::vector<FieldBase *> &fields = mesh_meta_data().get_fields();
  m_field_data_manager->reinitialize_removed_entity_field_data(rank, bucket_id, bucket_ord, fields);
}

void BulkData::destroy_bucket_callback(EntityRank rank, Bucket const& dying_bucket, unsigned capacity)
{
  // Remove destroyed bucket out of memoized get_buckets result, but
  // don't bother if the mesh is being destructed.
  const unsigned bucket_id = dying_bucket.bucket_id();

  if (!m_bucket_repository.being_destroyed()) {
    for (SelectorBucketMap::iterator itr = m_selector_to_buckets_map.begin(), end = m_selector_to_buckets_map.end();
         itr != end; ++itr) {
      Selector const& sel = itr->first.second;
      const EntityRank map_rank = itr->first.first;
      if (map_rank == rank && sel(dying_bucket)) {
#ifdef GATHER_GET_BUCKETS_METRICS
        ++m_num_buckets_removed_from_cache;
#endif
        TrackedBucketVector & cached_buckets = itr->second;
        TrackedBucketVector::iterator lb_itr = std::lower_bound(cached_buckets.begin(), cached_buckets.end(), bucket_id, BucketIdComparator());
        ThrowAssertMsg(lb_itr != cached_buckets.end() && (*lb_itr)->bucket_id() == bucket_id,
                       "Error, bucket id " << bucket_id << ":\n " << dying_bucket << "\nWas selected by selector " << sel << " but was not in the cache");
        cached_buckets.erase(lb_itr);
      }
    }
  }

  if (!m_keep_fields_updated) {
    return;
  }

  const std::vector<FieldBase*>&  fields = mesh_meta_data().get_fields();
  m_field_data_manager->deallocate_bucket_field_data(rank, bucket_id, capacity, fields);
}

void BulkData::update_field_data_states()
{
  const std::vector<FieldBase*> & field_set = mesh_meta_data().get_fields();

    for ( int i = 0 ; i < m_num_fields ; ) {
      const FieldBase & field = * field_set[i];
      const int outer_idx = i;
      const int num_state = field.number_of_states();
      i += num_state ;




      if (num_state > 1) {
        for ( int b = 0, be = field_set[outer_idx]->get_meta_data_for_field().size(); b < be; ++b) {
          if ( field_set[outer_idx]->get_meta_data_for_field()[b].m_bytes_per_entity > 0 ) {
            unsigned char* data_last = field_set[outer_idx]->get_meta_data_for_field()[b].m_data;
            for ( int s = 1; s < num_state; ++s ) {
              std::swap(field_set[outer_idx+s]->get_meta_data_for_field()[b].m_data, data_last);
            }
            field_set[outer_idx]->get_meta_data_for_field()[b].m_data = data_last;
          }
        }
      }

      for ( int s = 1; s < num_state; ++s )
      {
          m_field_data_manager->swap_fields(outer_idx+s, outer_idx);
      }

    }
}

void BulkData::reorder_buckets_callback(EntityRank rank, const std::vector<unsigned>& reorderedBucketIds)
{
  for (SelectorBucketMap::iterator itr = m_selector_to_buckets_map.begin(), end = m_selector_to_buckets_map.end();
       itr != end; ++itr) {
    TrackedBucketVector& cached_buckets = itr->second;
    std::sort(cached_buckets.begin(), cached_buckets.end(), BucketIdComparator());
  }

  if (!m_keep_fields_updated) {
    return;
  }

  const std::vector<FieldBase*>  fields = mesh_meta_data().get_fields();
  m_field_data_manager->reorder_bucket_field_data(rank, fields, reorderedBucketIds);
}

void BulkData::dump_all_mesh_info(std::ostream& out) const
{
  // Dump output for metadata first
  m_mesh_meta_data.dump_all_meta_info(out);

  out << "BulkData "
      << " info...\n";

  const FieldVector& all_fields = m_mesh_meta_data.get_fields();

  // Iterate all buckets for all ranks...
  const std::vector<std::string> & rank_names = m_mesh_meta_data.entity_rank_names();
  for (size_t i = 0, e = rank_names.size(); i < e; ++i) {
    EntityRank rank = static_cast<EntityRank>(i);
    out << "  All " << rank_names[i] << " entities:" << std::endl;

    const BucketVector& buckets = this->buckets(rank);
    BOOST_FOREACH(Bucket* bucket, buckets) {
      out << "    Found bucket " << bucket->bucket_id() << " with superset parts: { ";
      PartVector supersets;
      bucket->supersets(supersets);
      BOOST_FOREACH(Part* part, supersets) {
        out << part->name() << " ";
      }
      out << "}" << std::endl;

      for (size_t b_ord = 0, b_end = bucket->size(); b_ord < b_end; ++b_ord) {
        Entity entity = (*bucket)[b_ord];
        out << "      " << print_entity_key(m_mesh_meta_data, entity_key(entity)) << "(offset: " << entity.local_offset() << ")" << std::endl;

        // Print connectivity
        for (EntityRank r = stk::topology::NODE_RANK, re = static_cast<EntityRank>(rank_names.size()); r < re; ++r) {
          if (connectivity_map().valid(static_cast<EntityRank>(rank), r)) {
            out << "        Connectivity to " << rank_names[r] << std::endl;
            Entity const* entities = bucket->begin(b_ord, r);
            ConnectivityOrdinal const* ordinals = bucket->begin_ordinals(b_ord, r);
            const int num_conn         = bucket->num_connectivity(b_ord, r);
            for (int c_itr = 0; c_itr < num_conn; ++c_itr) {
              out << "          " << print_entity_key(m_mesh_meta_data, entity_key(entities[c_itr])) << "[" << ordinals[c_itr] << "]" << std::endl;
            }
          }
        }

        // Print field data
        if (m_num_fields > 0) {
          BOOST_FOREACH(FieldBase* field, all_fields) {

            if(static_cast<unsigned>(field->entity_rank()) != bucket->entity_rank()) continue;

            FieldMetaData field_meta_data = field->get_meta_data_for_field()[bucket->bucket_id()];

            unsigned data_size = field_meta_data.m_bytes_per_entity;
            if (data_size > 0) { // entity has this field?
              void* data = field_meta_data.m_data + field_meta_data.m_bytes_per_entity * b_ord;
              out << "        For field: " << *field << ", has data: ";
              field->print_data(out, data, data_size);
              out << std::endl;
            }
          }
        }
      }

    }
  }
}

void BulkData::set_relation_orientation(Entity from, Entity to, ConnectivityOrdinal to_ord, unsigned to_orientation)
{
  const EntityRank from_rank = entity_rank(from);
  const EntityRank to_rank   = entity_rank(to);

  Entity const*              fwd_rels  = begin(from, to_rank);
  ConnectivityOrdinal const* fwd_ords  = begin_ordinals(from, to_rank);
  Permutation *              fwd_perms = const_cast<Permutation*>(begin_permutations(from, to_rank));
  const int                  num_fwd   = num_connectivity(from, to_rank);

  Entity const*              back_rels  = begin(to, from_rank);
  ConnectivityOrdinal const* back_ords  = begin_ordinals(to, from_rank);
  Permutation *              back_perms = const_cast<Permutation*>(begin_permutations(to, from_rank));
  const int                  num_back   = num_connectivity(to,from_rank);

  // Find and change fwd connectivity
  for (int i = 0; i < num_fwd; ++i, ++fwd_rels, ++fwd_ords, ++fwd_perms) {
    // Allow clients to make changes to orientation
    // Orientations do not affect Relation ordering, so this is safe.
    if (*fwd_rels == to && *fwd_ords == to_ord) {
      *fwd_perms = static_cast<Permutation>(to_orientation);
    }
  }

  // Find and change back connectivity
  for (int i = 0; i < num_back; ++i, ++back_rels, ++back_ords, ++back_perms) {
    // Allow clients to make changes to orientation
    // Orientations do not affect Relation ordering, so this is safe.
    if (*back_rels == from && *back_ords == to_ord) {
      *back_perms = static_cast<Permutation>(to_orientation);
    }
  }
}

namespace {

void find_potential_upward_entities( const BulkData & mesh, Entity node, EntityRank rank, EntityVector & entities)
{
  entities.clear();

  // NOTE: it's adequate to just look at one node

  if ( mesh.connectivity_map().valid(stk::topology::NODE_RANK, rank) ) {
    entities.assign(mesh.begin(node, rank), mesh.end(node, rank));
  }
  else {
    // Have to go up to the elements, then back down to the target rank
    ThrowRequireMsg(mesh.connectivity_map().valid(stk::topology::NODE_RANK, stk::topology::ELEMENT_RANK), "ERROR, stk::mesh::get_connectivity was called but upward node->element connections are disabled. For get_connectivity to find non-stored connectivities, node->element connections must be enabled.");

    EntityVector elements;
    find_potential_upward_entities(mesh, node, stk::topology::ELEMENT_RANK, elements);

    for (unsigned i = 0, e = elements.size(); i < e; ++i) {
      entities.insert(entities.end(), mesh.begin(elements[i], rank), mesh.end(elements[i], rank));
    }

    std::sort(entities.begin(), entities.end());
    EntityVector::iterator new_end = std::unique(entities.begin(), entities.end());
    entities.erase(new_end, entities.end());
  }
}

template <bool GatherOrdinals, bool GatherPermutations>
size_t get_connectivity_impl(const BulkData & mesh,
                             Entity entity,
                             EntityRank to_rank,
                             EntityVector & entity_scratch_storage,
                             std::vector<ConnectivityOrdinal>* ordinal_scratch_storage = NULL,
                             std::vector<Permutation>* permutation_scratch_storage = NULL)
{
  ThrowAssert( !GatherOrdinals || ordinal_scratch_storage != NULL );
  ThrowAssert( !GatherPermutations || permutation_scratch_storage != NULL );

  const EntityRank source_rank = mesh.entity_rank(entity);

  // if the connectivity stored on the mesh, we shouldn't be calling this
  ThrowAssert( !mesh.connectivity_map().valid(source_rank, to_rank) );
  if ( source_rank >= to_rank || (mesh.mesh_meta_data().spatial_dimension() == 2 && to_rank == stk::topology::FACE_RANK) ) {
    // There is no possible connectivity
    return 0;
  }
  // go through the nodes, keeping in mind that we may be a node
  else {
    entity_scratch_storage.clear();
    if (GatherOrdinals)     { ordinal_scratch_storage->clear(); }
    if (GatherPermutations) { permutation_scratch_storage->clear(); }

    unsigned num_nodes;
    Entity const * nodes;
    if (source_rank == stk::topology::NODE_RANK) {
      num_nodes = 1;
      nodes     = &entity;
    }
    else {
      num_nodes = mesh.num_nodes(entity);
      nodes     = mesh.begin_nodes(entity);
    }

    if (num_nodes > 0) {
      EntityVector potential_upward_entities;
      find_potential_upward_entities( mesh, nodes[0], to_rank, potential_upward_entities);

      for (size_t i = 0u, e = potential_upward_entities.size(); i < e; ++i) {
        Entity potential_upward_entity = potential_upward_entities[i];
        Entity const * potential_sources = mesh.begin(potential_upward_entity, source_rank);
        ConnectivityOrdinal const * potential_ordinals = GatherOrdinals ? mesh.begin_ordinals(potential_upward_entity, source_rank) : NULL;
        Permutation const * potential_permutations = GatherPermutations ? mesh.begin_permutations(potential_upward_entity, source_rank) : NULL;
        const unsigned num_sources = mesh.num_connectivity(potential_upward_entity, source_rank);

        for (unsigned is=0u; is < num_sources; ++is) {
          if ( potential_sources[is] == entity) {
            entity_scratch_storage.push_back(potential_upward_entity);
            if (GatherOrdinals) { ordinal_scratch_storage->push_back(potential_ordinals[is]); }
            if (GatherPermutations && potential_permutations != NULL) { permutation_scratch_storage->push_back(potential_permutations[is]); }
          }
        }
      }
    }

    return entity_scratch_storage.size();
  }
}

static const bool GATHER_ORDINALS     = true;
static const bool GATHER_PERMUTATIONS = true;
static const bool DONT_GATHER_ORDINALS     = false;
static const bool DONT_GATHER_PERMUTATIONS = false;


} // namespace

size_t get_connectivity( const BulkData & mesh,
                         Entity entity,
                         EntityRank to_rank,
                         EntityVector & entity_scratch_storage)
{
  return get_connectivity_impl<DONT_GATHER_ORDINALS, DONT_GATHER_PERMUTATIONS>(mesh, entity, to_rank, entity_scratch_storage);
}

size_t get_connectivity( const BulkData & mesh,
                         Entity entity,
                         EntityRank to_rank,
                         EntityVector & entity_scratch_storage,
                         std::vector<ConnectivityOrdinal> & ordinals )
{
  return get_connectivity_impl<GATHER_ORDINALS, DONT_GATHER_PERMUTATIONS>(mesh, entity, to_rank, entity_scratch_storage, &ordinals);
}

size_t get_connectivity( const BulkData & mesh,
                         Entity entity,
                         EntityRank to_rank,
                         EntityVector & entity_scratch_storage,
                         std::vector<Permutation> & permutations )
{
  std::vector<ConnectivityOrdinal>* ignore = NULL;
  return get_connectivity_impl<DONT_GATHER_ORDINALS, GATHER_PERMUTATIONS>(mesh, entity, to_rank, entity_scratch_storage, ignore, &permutations);
}

size_t get_connectivity( const BulkData & mesh,
                         Entity entity,
                         EntityRank to_rank,
                         EntityVector & entity_scratch_storage,
                         std::vector<ConnectivityOrdinal> & ordinals,
                         std::vector<Permutation> & permutations )
{
  return get_connectivity_impl<GATHER_ORDINALS, GATHER_PERMUTATIONS>(mesh, entity, to_rank, entity_scratch_storage, &ordinals, &permutations);
}


void BulkData::reserve_relation(Entity entity, const unsigned num)
{
  if (num == 0 && aux_relations(entity).empty()) {
    RelationVector tmp;
    aux_relations(entity).swap(tmp); // clear memory of m_relations.
  }
  else {
    aux_relations(entity).reserve(num);
  }
}

void BulkData::erase_and_clear_if_empty(Entity entity, RelationIterator rel_itr)
{
  ThrowAssert(!impl::internal_is_handled_generically(rel_itr->getRelationType()));

  RelationVector& aux_rels = aux_relations(entity);
  aux_rels.erase(aux_rels.begin() + (rel_itr - aux_rels.begin())); // Need to convert to non-const iterator

  if (aux_rels.empty()) {
    reserve_relation(entity, 0);
  }
}

void BulkData::internal_verify_initialization_invariant(Entity entity)
{
#ifndef NDEBUG
  int my_global_id = global_id(entity);
  EntityKey my_key = entity_key(entity);
#endif
  ThrowAssert ( !(my_global_id < 0 && my_key.id() == static_cast<EntityId>(my_global_id)) &&
                !(my_global_id > 0 && my_key.id() != static_cast<EntityId>(my_global_id)) );
}

BucketVector const& BulkData::get_buckets(EntityRank rank, Selector const& selector) const
{
#ifdef GATHER_GET_BUCKETS_METRICS
  ++m_num_memoized_get_buckets_calls;
#endif

  std::pair<EntityRank, Selector> search_item = std::make_pair(rank, selector);
  SelectorBucketMap::iterator fitr =
    m_selector_to_buckets_map.find(search_item);
  if (fitr != m_selector_to_buckets_map.end()) {
#ifdef GATHER_GET_BUCKETS_METRICS
    std::pair<size_t, size_t> & data = m_selector_to_count_map[search_item];
    ++data.first;
    data.second += buckets(rank).size();
#endif
    TrackedBucketVector const& rv = fitr->second;
    return reinterpret_cast<BucketVector const&>(rv);
  }
  else {
    BucketVector const& all_buckets_for_rank = buckets(rank); // lots of potential side effects! Need to happen before map insertion
    std::pair<SelectorBucketMap::iterator, bool> insert_rv =
      m_selector_to_buckets_map.insert(std::make_pair( std::make_pair(rank, selector), TrackedBucketVector() ));
    ThrowAssertMsg(insert_rv.second, "Should not have already been in map");
    TrackedBucketVector& map_buckets = insert_rv.first->second;
    for (size_t i = 0, e = all_buckets_for_rank.size(); i < e; ++i) {
      if (selector(*all_buckets_for_rank[i])) {
        map_buckets.push_back(all_buckets_for_rank[i]);
      }
    }

#ifdef GATHER_GET_BUCKETS_METRICS
    m_selector_to_count_map[search_item] = std::make_pair(0, 0); // (hits, bucket traversals saved)
#endif

    return reinterpret_cast<BucketVector const&>(map_buckets);
  }
}

void BulkData::get_buckets(EntityRank rank, Selector const& selector, BucketVector & output_buckets) const
{
#ifdef GATHER_GET_BUCKETS_METRICS
  ++m_num_non_memoized_get_buckets_calls;
#endif

  output_buckets.clear();

  BucketVector const& all_buckets_for_rank = buckets(rank);
  for (size_t i = 0, e = all_buckets_for_rank.size(); i < e; ++i) {
    if (selector(*all_buckets_for_rank[i])) {
      output_buckets.push_back(all_buckets_for_rank[i]);
    }
  }
}

void BulkData::get_entities(EntityRank rank, Selector const& selector, EntityVector& output_entities) const {
  output_entities.clear();
  const stk::mesh::BucketVector &bucket_ptrs = get_buckets(rank, selector);
  for(size_t ib=0, ib_end=bucket_ptrs.size(); ib<ib_end; ++ib) {
     const stk::mesh::Bucket *bucket = bucket_ptrs[ib];
     for(size_t iobj=0, iobj_end=bucket->size(); iobj<iobj_end; ++iobj) {
       output_entities.push_back((*bucket)[iobj]);
     }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

//----------------------------------------------------------------------

void BulkData::require_valid_relation( const char action[] ,
                                       const BulkData & mesh ,
                                       const Entity e_from ,
                                       const Entity e_to )
{
  const bool error_type      = mesh.entity_rank(e_from) <= mesh.entity_rank(e_to);
  const bool error_nil_from  = !mesh.is_valid(e_from);
  const bool error_nil_to    = !mesh.is_valid(e_to);

  if ( error_type || error_nil_from || error_nil_to ) {
    std::ostringstream msg ;

    msg << "Could not " << action << " relation from entity "
        << print_entity_key(MetaData::get(mesh), mesh.entity_key(e_from)) << " to entity "
        << print_entity_key(MetaData::get(mesh), mesh.entity_key(e_to)) << "\n";

    ThrowErrorMsgIf( error_nil_from  || error_nil_to,
                     msg.str() << ", entity was destroyed");
    ThrowErrorMsgIf( error_type, msg.str() <<
                     "A relation must be from higher to lower ranking entity");
  }
}


//----------------------------------------------------------------------
bool BulkData::internal_declare_relation(Entity e_from, Entity e_to,
                                         RelationIdentifier local_id,
                                         unsigned sync_count, bool is_back_relation,
                                         Permutation permut)
{
    INCREMENT_MODIFICATION_COUNTER(INTERNAL, DECLARE_RELATION);

  TraceIfWatching("stk::mesh::BuilkData::internal_declare_relation", LOG_ENTITY, entity_key(e_from));

  const MeshIndex& idx = mesh_index(e_from);

  bool modified = idx.bucket->declare_relation(idx.bucket_ordinal, e_to, static_cast<ConnectivityOrdinal>(local_id),
                                               permut);

  if (modified) {
    set_synchronized_count( e_from, sync_count );
  }
  return modified;
}

void BulkData::declare_relation( Entity e_from ,
                                 Entity e_to ,
                                 const RelationIdentifier local_id ,
                                 Permutation permut)
{
  INCREMENT_MODIFICATION_COUNTER(PUBLIC, DECLARE_RELATION);
  OrdinalVector ordinal_scratch;
  PartVector part_scratch;
  internal_declare_relation(e_from, e_to, local_id, permut, ordinal_scratch, part_scratch);
}

void BulkData::declare_relation( Entity e_from ,
                                 Entity e_to ,
                                 const RelationIdentifier local_id ,
                                 Permutation permut,
                                 OrdinalVector& ordinal_scratch,
                                 PartVector& part_scratch)
{
    INCREMENT_MODIFICATION_COUNTER(PUBLIC, DECLARE_RELATION);
    internal_declare_relation(e_from, e_to, local_id, permut, ordinal_scratch, part_scratch);
}

void BulkData::internal_declare_relation( Entity e_from ,
                                 Entity e_to ,
                                 const RelationIdentifier local_id ,
                                 Permutation permut,
                                 OrdinalVector& ordinal_scratch,
                                 PartVector& part_scratch)
{
  TraceIfWatching("stk::mesh::BulkData::declare_relation", LOG_ENTITY, entity_key(e_from));
  TraceIfWatchingDec("stk::mesh::BulkData::declare_relation", LOG_ENTITY, entity_key(e_to), 1);
  DiagIfWatching(LOG_ENTITY, entity_key(e_from),
                 "from: " << entity_key(e_from) << ";  " <<
                 "to: " << entity_key(e_to) << ";  " <<
                 "id: " << local_id);
  DiagIfWatching(LOG_ENTITY, entity_key(e_to),
                 "from: " << entity_key(e_from) << ";  " <<
                 "to: " << entity_key(e_to) << ";  " <<
                 "id: " << local_id);

  require_ok_to_modify();

  require_valid_relation( "declare" , *this , e_from , e_to );

  // TODO: Don't throw if exact relation already exists, that should be a no-op.
  // Should be an exact match if relation of local_id already exists (e_to should be the same).
  bool is_converse = false;
  bool caused_change_fwd = internal_declare_relation(e_from, e_to, local_id, m_sync_count,
                                                     is_converse, permut);

  //TODO: check connectivity map
  // Relationships should always be symmetrical
  if ( caused_change_fwd ) {

    const bool higher_order_relation = stk::topology::ELEMENT_RANK < entity_rank(e_from);
    if (    higher_order_relation
         || m_bucket_repository.connectivity_map()(entity_rank(e_to), entity_rank(e_from)) != stk::mesh::INVALID_CONNECTIVITY_TYPE
       )
    {
      // the setup for the converse relationship works slightly differently
      is_converse = true;
      internal_declare_relation(e_to, e_from, local_id, m_sync_count, is_converse, permut );
    }
  }

  // It is critical that the modification be done AFTER the relations are
  // added so that the propagation can happen correctly.
  if ( caused_change_fwd ) {
    this->mark_entity_and_upward_related_entities_as_modified(e_to);
    this->mark_entity_and_upward_related_entities_as_modified(e_from);
  }

  OrdinalVector empty ;

  // Deduce and set new part memberships:
  ordinal_scratch.clear();

  induced_part_membership(*this, e_from, empty, entity_rank(e_to), ordinal_scratch );

  PartVector emptyParts;
  part_scratch.clear();
  for(unsigned ipart=0; ipart<ordinal_scratch.size(); ++ipart) {
    part_scratch.push_back(&m_mesh_meta_data.get_part(ordinal_scratch[ipart]));
  }

  internal_change_entity_parts( e_to , part_scratch , emptyParts );
}

//----------------------------------------------------------------------

void BulkData::internal_declare_relation( Entity entity ,
                                 const std::vector<Relation> & rel )
{
  OrdinalVector ordinal_scratch;
  PartVector part_scratch;
  internal_declare_relation(entity, rel, ordinal_scratch, part_scratch);
}

void BulkData::internal_declare_relation( Entity entity ,
                                 const std::vector<Relation> & rel,
                                 OrdinalVector& ordinal_scratch,
                                 PartVector& part_scratch )
{
  require_ok_to_modify();

  const unsigned erank = entity_rank(entity);

  std::vector<Relation>::const_iterator i ;
  for ( i = rel.begin() ; i != rel.end() ; ++i ) {
    Entity e = i->entity();
    const unsigned n = i->relation_ordinal();
    const Permutation permut = static_cast<Permutation>(i->getOrientation());
    if ( entity_rank(e) < erank ) {
      internal_declare_relation( entity , e , n, permut, ordinal_scratch, part_scratch );
    }
    else if ( erank < entity_rank(e) ) {
      internal_declare_relation( e , entity , n, permut, ordinal_scratch, part_scratch );
    }
    else {
      ThrowErrorMsg("Given entities of the same entity rank. entity is " <<
                    identifier(entity));
    }
  }
}

//----------------------------------------------------------------------

bool BulkData::destroy_relation( Entity e_from ,
                                 Entity e_to,
                                 const RelationIdentifier local_id )
{
    INCREMENT_MODIFICATION_COUNTER(PUBLIC, DESTROY_RELATION);
    return internal_destroy_relation(e_from, e_to,  local_id);
}

bool BulkData::internal_destroy_relation( Entity e_from ,
                                 Entity e_to,
                                 const RelationIdentifier local_id )
{
    INCREMENT_MODIFICATION_COUNTER(INTERNAL, DESTROY_RELATION);

  TraceIfWatching("stk::mesh::BulkData::destroy_relation", LOG_ENTITY, entity_key(e_from));
  TraceIfWatchingDec("stk::mesh::BulkData::destroy_relation", LOG_ENTITY, entity_key(e_to), 1);
  DiagIfWatching(LOG_ENTITY, entity_key(e_from),
                 "from: " << entity_key(e_from) << ";  " <<
                 "to: " << entity_key(e_to) << ";  " <<
                 "id: " << local_id);
  DiagIfWatching(LOG_ENTITY, entity_key(e_to),
                 "from: " << entity_key(e_from) << ";  " <<
                 "to: " << entity_key(e_to) << ";  " <<
                 "id: " << local_id);

  require_ok_to_modify();

  require_valid_relation( "destroy" , *this , e_from , e_to );

  const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());
  const EntityRank e_to_entity_rank = entity_rank(e_to);
  const EntityRank e_from_entity_rank = entity_rank(e_from);

  //------------------------------
  // When removing a relationship may need to
  // remove part membership and set field relation pointer to NULL

  m_check_invalid_rels = false; // OK to have gaps when deleting

  if ( parallel_size() < 2 || entity_comm_map_shared(entity_key(e_to)).empty() ) {

    //------------------------------
    // 'keep' contains the parts deduced from kept relations
    // 'del'  contains the parts deduced from deleted relations
    //        that are not in 'keep'
    // Only remove these part memberships the entity is not shared.
    // If the entity is shared then wait until modificaton_end_synchronize.
    //------------------------------

    OrdinalVector del, keep, empty;

    // For all relations that are *not* being deleted, add induced parts for
    // these relations to the 'keep' vector
    {
      EntityVector temp_entities;
      std::vector<ConnectivityOrdinal> temp_ordinals;
      Entity const* rel_entities = NULL;
      ConnectivityOrdinal const* rel_ordinals = NULL;
      int num_rels = 0;
      for (EntityRank irank = static_cast<EntityRank>(e_to_entity_rank + 1); irank < end_rank; ++irank) {
        if (connectivity_map().valid(e_to_entity_rank, irank)) {
          num_rels     = num_connectivity(e_to, irank);
          rel_entities = begin(e_to, irank);
          rel_ordinals = begin_ordinals(e_to, irank);
        }
        else {
          num_rels     = get_connectivity(*this, e_to, irank, temp_entities, temp_ordinals);
          rel_entities = &*temp_entities.begin();
          rel_ordinals = &*temp_ordinals.begin();
        }

        for (int j = 0; j < num_rels; ++j) {
          ThrowAssertMsg(is_valid(rel_entities[j]), "Error, entity " << e_to.local_offset() << " has invalid back-relation for ordinal: "
                         << rel_ordinals[j] << " to rank: " << irank << ", target entity is: " << rel_entities[j].local_offset());
          if ( !(rel_entities[j] == e_from && rel_ordinals[j] == static_cast<ConnectivityOrdinal>(local_id) ) )
          {
            induced_part_membership(*this, rel_entities[j], empty, e_to_entity_rank, keep,
                                    false /*Do not look at supersets*/);
          }
        }
      }
    }

    // Find the relation this is being deleted and add the parts that are
    // induced from that relation (and that are not in 'keep') to 'del'
    {
      size_t num_rels = num_connectivity(e_from, e_to_entity_rank);
      Entity const *rel_entities = begin(e_from, e_to_entity_rank);
      ConnectivityOrdinal const *rel_ordinals = begin_ordinals(e_from, e_to_entity_rank);
      for (size_t j = 0; j < num_rels; ++j)
      {
        if ( rel_entities[j] == e_to && rel_ordinals[j] == static_cast<ConnectivityOrdinal>(local_id) )
        {
          induced_part_membership(*this, e_from, keep, e_to_entity_rank, del,
                                  false /*Do not look at supersets*/);
          break; // at most 1 relation can match our specification
        }
      }
    }

    if ( !del.empty() ) {

      PartVector delParts, emptyParts;
      delParts.reserve(del.size());
      for(unsigned ipart=0; ipart<del.size(); ++ipart) {
        delParts.push_back(&m_mesh_meta_data.get_part(del[ipart]));
      }

      internal_change_entity_parts( e_to , emptyParts , delParts );
    }
  }

  //delete relations from the entities
  bool caused_change_fwd = bucket(e_from).destroy_relation(e_from, e_to, local_id);

  // Relationships should always be symmetrical
  if ( caused_change_fwd &&
       (e_to_entity_rank > stk::topology::ELEMENT_RANK || e_from_entity_rank > stk::topology::ELEMENT_RANK ||
        connectivity_map().valid(entity_rank(e_to), entity_rank(e_from))) ) {
    bucket(e_to).destroy_relation(e_to, e_from, local_id);
  }

  // It is critical that the modification be done AFTER the relations are
  // changed so that the propagation can happen correctly.
  if ( caused_change_fwd ) {
    this->mark_entity_and_upward_related_entities_as_modified(e_to);
    this->mark_entity_and_upward_related_entities_as_modified(e_from);
  }

  m_check_invalid_rels = true;

  return caused_change_fwd;
}



//----------------------------------------------------------------------
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//ParallelVerify
//----------------------------------------------------------------------

namespace {

//#define DEBUG_PRINT_COMM_LIST
//#define DEBUG_PRINT_COMM_LIST_UNPACK

#ifdef DEBUG_PRINT_COMM_LIST

// Very, very handy for debugging parallel resolution...
static int s_step = 0;

void par_verify_print_comm_list( const BulkData & mesh , bool doit, const std::string message )
{
  ++s_step;
  if ( doit ) {
    std::ostringstream file ;
    file << "comm-list." << s_step << "." << mesh.parallel_rank() << ".dat";
    std::ofstream fout(file.str().c_str());
    std::ostringstream msg ;
    msg << message;
    msg << std::endl ;

    for ( EntityCommListInfoVector::const_iterator
          i =  mesh.comm_list().begin() ;
          i != mesh.comm_list().end() ; ++i ) {

      Entity entity = i->entity;
      msg << "S< " << s_step << " > P" << mesh.parallel_rank() << ": " ;

      print_entity_key( msg , MetaData::get(mesh) , i->key );

      msg << " owner(" << i->owner << ")" ;

      if ( !mesh.is_valid(entity) ) { msg << " del" ; }
      else if ( Modified == mesh.state(entity) ) { msg << " mod" ; }
      else { msg << "    " ; }

      for ( PairIterEntityComm ec = mesh.entity_comm_map(i->key); ! ec.empty() ; ++ec ) {
        msg << " (" << ec->ghost_id << "," << ec->proc << ")" ;
      }
      msg << std::endl ;
    }

    fout << msg.str();
  }
}

#endif

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::change_entity_owner( const std::vector<EntityProc> & arg_change,
                                    bool regenerate_aura,
                                    modification_optimization mod_optimization )
{
    INCREMENT_MODIFICATION_COUNTER(PUBLIC, CHANGE_ENTITY_OWNER);
    internal_change_entity_owner(arg_change, regenerate_aura, mod_optimization);
    internal_modification_end_for_change_entity_owner(regenerate_aura, mod_optimization);
}

void BulkData::internal_compute_proposed_owned_closure_count(const std::vector<EntityProc> & local_change,
                                                             const std::vector<EntityProc> & shared_change,
                                                             const std::vector<EntityProc> & ghosted_change,
                                                             std::vector<uint16_t> & new_closure_count) const
{
    new_closure_count.assign(m_closure_count.begin(),m_closure_count.end());
    for (size_t i=0 ; i<local_change.size() ; ++i) {
        Entity locally_owned_entity = local_change[i].first;
        // Decrement because we're giving ownership away.
        --new_closure_count[locally_owned_entity.local_offset()];

        // Decrement all downwardly connected entities because this locally_owned_entity ownership is being given away and therefore it will not be locally owned after change_entity_owner.
        for (EntityRank rank = stk::topology::NODE_RANK, end_rank = entity_rank(locally_owned_entity); rank < end_rank; ++rank) {
            unsigned num = num_connectivity(locally_owned_entity,rank);
            Entity const * entities = begin(locally_owned_entity,rank);
            for (unsigned j =0; j<num; ++j) {
                --new_closure_count[entities[j].local_offset()];
            }
        }
    }
    for (size_t i=0 ; i<shared_change.size() ; ++i) {
        Entity remotely_owned_entity = shared_change[i].first;
        int new_owner_proc = shared_change[i].second;
        if (new_owner_proc == parallel_rank()) {
            // Increment because we're being given ownership.
            ++new_closure_count[remotely_owned_entity.local_offset()];

            // Increment all downwardly connected entities because this remotely_owned_entity ownership is being given to us and therefore it will be locally owned after change_entity_owner.
            for (EntityRank rank = stk::topology::NODE_RANK, end_rank = entity_rank(remotely_owned_entity); rank < end_rank; ++rank) {
                unsigned num = num_connectivity(remotely_owned_entity,rank);
                Entity const * entities = begin(remotely_owned_entity,rank);
                for (unsigned j =0; j<num; ++j) {
                    ++new_closure_count[entities[j].local_offset()];
                }
            }
        }
    }
    for (size_t i=0 ; i<ghosted_change.size() ; ++i) {
        Entity remotely_owned_entity = ghosted_change[i].first;
        int new_owner_proc = ghosted_change[i].second;
        if (new_owner_proc == parallel_rank()) {
            // Increment because we're being given ownership.
            ++new_closure_count[remotely_owned_entity.local_offset()];

            // Increment all downwardly connected entities because this remotely_owned_entity ownership is being given to us and therefore it will be locally owned after change_entity_owner.
            for (EntityRank rank = stk::topology::NODE_RANK, end_rank = entity_rank(remotely_owned_entity); rank < end_rank; ++rank) {
                unsigned num = num_connectivity(remotely_owned_entity,rank);
                Entity const * entities = begin(remotely_owned_entity,rank);
                for (unsigned j =0; j<num; ++j) {
                    ++new_closure_count[entities[j].local_offset()];
                }
            }
        }
    }
}

namespace {
void print_entity_to_dependent_processors_map(std::string title, int parallel_rank, const BulkData::NodeToDependentProcessorsMap & entity_to_dependent_processors_map)
{
#ifndef NDEBUG
  {
      std::ostringstream oss;
      oss << "\n";
      BulkData::NodeToDependentProcessorsMap::const_iterator dep_proc_it = entity_to_dependent_processors_map.begin();
      for (; dep_proc_it != entity_to_dependent_processors_map.end() ; ++dep_proc_it ) {
          EntityKey key = dep_proc_it->first;
          const std::set<int> & sharing_procs = dep_proc_it->second;
          oss << title << " P" << parallel_rank << ": entity_to_dependent_processors_map[" << key << "] = { ";
          for (std::set<int>::const_iterator it=sharing_procs.begin() ; it!=sharing_procs.end() ; ++it) {
              oss << *it << " ";
          }
          oss << "}" << std::endl;
          std::cout << oss.str() << std::flush;
      }
  }
#endif // NDEBUG

}

} // namespace

void BulkData::internal_communicate_entity_to_dependent_processors_map(
        NodeToDependentProcessorsMap & entity_to_dependent_processors_map)
{
  print_entity_to_dependent_processors_map("BEFORE",parallel_rank(),entity_to_dependent_processors_map);
  // Now we have updated sharing lists for every soon-to-be not-owned node.
  // Lets talk with our neighbors!
  const int p_size = parallel_size();
  CommAll comm( parallel() );
  std::vector<int> procs ;

  // pack and communicate sharing information to all
  // processes that need to know about the entity
  for ( int phase = 0; phase < 2; ++phase) {
    NodeToDependentProcessorsMap::const_iterator map_it = entity_to_dependent_processors_map.begin();
    for ( ; map_it != entity_to_dependent_processors_map.end() ; ++map_it ) {
        EntityKey node_key = map_it->first;
        const std::set<int> & sharing_procs = map_it->second;
        std::set<int>::const_iterator proc_it = sharing_procs.begin();
        for ( ; proc_it != sharing_procs.end() ; ++proc_it) {
            if (*proc_it == parallel_rank()) { continue; }
            comm.send_buffer(*proc_it)
                                    .pack<EntityKey>( node_key )
                                    .pack<int>( sharing_procs.size() );
            std::set<int>::const_iterator proc_it_inner = sharing_procs.begin();
            for ( ; proc_it_inner != sharing_procs.end() ; ++proc_it_inner) {
                comm.send_buffer(*proc_it).pack<int>( *proc_it_inner );
            }
        }
    }
    if (phase == 0) { // allocation phase
      comm.allocate_buffers( p_size / 4 , 0 );
    }
    else { // communication phase
      comm.communicate();
    }
  }

  // unpack communicated owner information into the
  // ghosted and shared change vectors.
  for ( int ip = 0 ; ip < p_size ; ++ip ) {
    CommBuffer & buf = comm.recv_buffer( ip );
    while ( buf.remaining() ) {
      EntityKey key ;
      int num_sharing_procs = -1;
      buf.unpack<EntityKey>( key )
         .unpack<int>( num_sharing_procs );
      for (int i=0 ; i < num_sharing_procs ; ++i) {
          int sharing_proc = -1;
          buf.unpack<int>( sharing_proc );
          entity_to_dependent_processors_map[key].insert(sharing_proc);
      }
    }
  }
  print_entity_to_dependent_processors_map("AFTER",parallel_rank(),entity_to_dependent_processors_map);

}

void internal_print_comm_list(std::string title, BulkData & mesh)
{
#ifndef NDEBUG
  std::ostringstream oss;
  oss << "\n";
  const EntityCommListInfoVector & entity_comm_list = mesh.comm_list();
  int index=0;
  for(EntityCommListInfoVector::const_iterator i = entity_comm_list.begin(); i != entity_comm_list.end(); ++i) {
      oss << title << " P" << mesh.parallel_rank() << ": m_entity_comm_list[" << index << "] owner=" << i->owner << ", key=" << i->key << ", owner_rank=" << i->entity_comm->owner_rank << ", comm_map = { ";
      const EntityCommInfoVector & comm_map = i->entity_comm->comm_map;
      for (EntityCommInfoVector::const_iterator comm_map_it = comm_map.begin() ; comm_map_it != comm_map.end() ; ++comm_map_it ) {
          unsigned ghost_id = comm_map_it->ghost_id;
          int proc = comm_map_it->proc;
          oss << "(" << (ghost_id == 0 ? "SHARED" : "AURA") << ", " << proc << ") " ;
      }
      oss << " }\n";
      ++index;
  }
  std::cout << oss.str() << std::flush;
#endif // NDEBUG
}

void BulkData::internal_calculate_sharing(const std::vector<EntityProc> & local_change,
                         const std::vector<EntityProc> & shared_change,
                         const std::vector<EntityProc> & ghosted_change,
                         NodeToDependentProcessorsMap & entity_to_dependent_processors_map)

{
  internal_print_comm_map("BEFORE");
  internal_print_comm_list("BEFORE", *this);
  typedef std::set<EntityProc,EntityLess> EntityProcSet;

  entity_to_dependent_processors_map.clear();
  for ( std::vector<EntityProc>::const_iterator
        owned_giving_away_eproc = local_change.begin() ; owned_giving_away_eproc != local_change.end() ; ++owned_giving_away_eproc ) {
      stk::mesh::impl::internal_get_processor_dependencies_shared_or_ghosted(*this, *owned_giving_away_eproc , entity_to_dependent_processors_map);
  }
  for (std::vector<EntityProc>::const_iterator
          shared_giving_away_eproc = shared_change.begin() ; shared_giving_away_eproc != shared_change.end() ; ++shared_giving_away_eproc ) {
      stk::mesh::impl::internal_get_processor_dependencies_shared_or_ghosted(*this, *shared_giving_away_eproc, entity_to_dependent_processors_map );
  }
  for (std::vector<EntityProc>::const_iterator
          ghosted_giving_away_eproc = ghosted_change.begin() ; ghosted_giving_away_eproc != ghosted_change.end() ; ++ghosted_giving_away_eproc ) {
      stk::mesh::impl::internal_get_processor_dependencies_shared_or_ghosted(*this, *ghosted_giving_away_eproc, entity_to_dependent_processors_map );
  }

#ifndef NDEBUG
  for (size_t i=1 ; i<m_closure_count.size() ; ++i) {
      Entity entity;
      entity.set_local_offset(i);
      if (!is_valid(entity)) { continue; }
      EntityKey key = entity_key(entity);
      std::ostringstream oss;
      oss << "P" << parallel_rank() << ": m_closure_count[" << i << "] " << key << " = " << m_closure_count[i] << std::endl;
      std::cout << oss.str() << std::flush;
  }
#endif // NDEBUG
  std::vector<uint16_t> proposed_closure_count;
  internal_compute_proposed_owned_closure_count(local_change, shared_change, ghosted_change, proposed_closure_count);
#ifndef NDEBUG
  for (size_t i=1 ; i<proposed_closure_count.size() ; ++i) {
      Entity entity;
      entity.set_local_offset(i);
      if (!is_valid(entity)) { continue; }
      EntityKey key = entity_key(entity);
      std::ostringstream oss;
      oss << "P" << parallel_rank() << ": proposed_closure_count[" << i << "] " << key << " = " << proposed_closure_count[i] << std::endl;
      std::cout << oss.str() << std::flush;
  }
#endif // NDEBUG

  update_dependent_processor_map_from_closure_count(proposed_closure_count, entity_to_dependent_processors_map);

  internal_communicate_entity_to_dependent_processors_map(entity_to_dependent_processors_map);
}

void BulkData::update_dependent_processor_map_from_closure_count(const std::vector<uint16_t> & proposed_closure_count,
                                                                 NodeToDependentProcessorsMap & entity_to_dependent_processors_map)
{
  // Do I depend on any of these entities after all the ownership changes?
  NodeToDependentProcessorsMap::iterator map_it = entity_to_dependent_processors_map.begin();
  for (; map_it != entity_to_dependent_processors_map.end() ; ++map_it) {
      Entity entity_in_send_closure = get_entity(map_it->first);
      if (proposed_closure_count[entity_in_send_closure.local_offset()] > static_cast<uint16_t>(0u)) {
          map_it->second.insert(parallel_rank());
          if (proposed_closure_count[entity_in_send_closure.local_offset()] == static_cast<uint16_t>(1u)) {
              std::ostringstream oss;
              oss << "P" << parallel_rank() << ": Warning: orphaned node = " << entity_in_send_closure << " is locally owned but no owned entities have downward relations to it on this processor!  This will likely brake I/O." << std::endl;
              std::cout << oss.str() << std::flush;
          }
      }
      else { // 0u
          map_it->second.erase(parallel_rank());
          // This entity can be removed from globally_shared part on my processor
      }
  }
}

void BulkData::internal_apply_node_sharing(const NodeToDependentProcessorsMap & entity_sharing_map)
{
    std::vector<stk::mesh::Entity> modified_entities;
    // update sharing!
    NodeToDependentProcessorsMap::const_iterator map_it = entity_sharing_map.begin();
    for ( ; map_it != entity_sharing_map.end() ; ++map_it ) {
        EntityKey key = map_it->first;
        entity_comm_map_erase(  key, shared_ghosting() );
        Entity entity = get_entity(key);
        if (!is_valid(entity)) { continue; }
        const std::set<int> & sharing_procs = map_it->second;
        if (sharing_procs.find(parallel_rank()) == sharing_procs.end()) {
            // Skip adding sharing if my processor is not listed in sharing_procs.
            continue;
        }
        ThrowAssert(entity_rank(entity) != stk::topology::ELEMENT_RANK); // TODO:  remove
        std::set<int>::const_iterator set_it = sharing_procs.begin();
        for (; set_it != sharing_procs.end() ; ++set_it ) {
            if (*set_it == parallel_rank()) { continue; }
            mark_entity(entity, IS_SHARED);
            m_add_node_sharing_called = true;
            // Cannot have both sharing and ghosting of a node to same proc.
            entity_comm_map_erase(key, EntityCommInfo(stk::mesh::BulkData::AURA, *set_it));
            entity_comm_map_insert(entity, EntityCommInfo(stk::mesh::BulkData::SHARED, *set_it));
            modified_entities.push_back(entity);
        }
        if (!bucket(entity).owned()) { continue; }
        PartVector globally_shared_pv;
        globally_shared_pv.push_back(&m_mesh_meta_data.globally_shared_part());
        PartVector empty_pv;
        const bool more_than_owner_in_sharing_list = (sharing_procs.size() > 1);
        if (more_than_owner_in_sharing_list) {
            internal_verify_and_change_entity_parts(entity, globally_shared_pv, empty_pv);
        }
        else {
            internal_verify_and_change_entity_parts(entity, empty_pv, globally_shared_pv);
        }
    }
    update_comm_list_based_on_changes_in_comm_map();

    update_comm_list(modified_entities);

    internal_print_comm_map("AFTER");
    internal_print_comm_list("AFTER", *this);
}

void BulkData::internal_print_comm_map (std::string title)
{
#ifndef NDEBUG
    std::ostringstream oss;
    for (EntityRank rank = stk::topology::NODE_RANK; rank != stk::topology::END_RANK; rank++){
        const_entity_iterator entity_it = begin_entities(rank);
        for(; entity_it != end_entities(rank); entity_it ++){
            EntityKey key = entity_it->first;
            PairIterEntityComm comm_pit = entity_comm_map(key);
            oss << title << " P" << parallel_rank() << ": m_entity_comm_map[" << key << "] = { ";
            for(; comm_pit.first != comm_pit.second; ++comm_pit){
                unsigned ghost_id = comm_pit.first->ghost_id;
                int proc = comm_pit.first->proc;
                oss << "(" << (ghost_id == 0 ? "SHARED" : "AURA") << ", " << proc << ") " ;
            }
            oss << "}" << std::endl;
        }
    }
    std::cout << oss.str() << std::flush;
#endif //NDEBUG
}

void BulkData::internal_change_entity_owner( const std::vector<EntityProc> & arg_change,
                                             bool regenerate_aura,
                                             modification_optimization mod_optimization )
{
  INCREMENT_MODIFICATION_COUNTER(INTERNAL, CHANGE_ENTITY_OWNER);

  Trace_("stk::mesh::BulkData::change_entity_owner");
  DiagIf(LOG_ENTITY, "arg_change: " << arg_change);

  const bool modStatus = modification_begin("change_entity_owner");
  ThrowRequireMsg(modStatus, "BulkData::change_entity_owner() must not be called from within a modification cycle.");

  require_ok_to_modify();

  const MetaData  & meta = m_mesh_meta_data ;
  const int       p_rank = m_parallel_rank ;
  const int       p_size = m_parallel_size ;
  ParallelMachine p_comm = m_parallel_machine ;

  //------------------------------
  // Verify the input changes, generate a clean local change list, and
  // generate the remote change list so that all processes know about
  // pending changes.

  std::vector<EntityProc> local_change( arg_change );

  // Parallel synchronous clean up and verify the requested changes:
  impl::internal_clean_and_verify_parallel_change( *this , local_change );

  //----------------------------------------
  // Parallel synchronous determination of changing shared and ghosted.

  // The two vectors below will contain changes to ghosted and shared
  // entities on this process coming from change-entity-owner requests
  // on other processes.
  std::vector<EntityProc> ghosted_change ;
  std::vector<EntityProc> shared_change ;

  impl::internal_generate_parallel_change_lists( *this , local_change ,
                            shared_change , ghosted_change );

#undef CHANGE_ENTITY_OWNER_UPDATES_COMM_INFO
#ifdef CHANGE_ENTITY_OWNER_UPDATES_COMM_INFO
  NodeToDependentProcessorsMap node_to_dependent_processors_map;
  internal_calculate_sharing(local_change, shared_change, ghosted_change, node_to_dependent_processors_map);
#endif
  //------------------------------
  // Have enough information to delete all effected ghosts.
  // If the closure of a ghost contains a changing entity
  // then that ghost must be deleted.
  // Request that all ghost entities in the closure of the ghost be deleted.

  StoreEntityProcInSet store_entity_proc_in_set(*this);

  // Compute the closure of all the locally changing entities
  for ( std::vector<EntityProc>::iterator
        i = local_change.begin() ; i != local_change.end() ; ++i ) {
      store_entity_proc_in_set.target = i->second;
      impl::VisitClosureGeneral(*this,i->first,store_entity_proc_in_set,store_entity_proc_in_set);
  }
  std::set<EntityProc,EntityLess> & send_closure = store_entity_proc_in_set.entity_proc_set;


  // Calculate all the ghosts that are impacted by the set of ownership
  // changes. We look at ghosted, shared, and local changes looking for ghosts
  // that are either in the closure of the changing entity, or have the
  // changing entity in their closure. All modified ghosts will be removed.
  {
      impl::OnlyVisitGhostsOnce only_visit_ghosts_once(*this);
      StoreEntityKeyInSet store_entity_key(*this);
      for ( std::vector<EntityProc>::const_iterator i = ghosted_change.begin() ; i != ghosted_change.end() ; ++i) {
          impl::VisitAuraClosureGeneral(*this,i->first,store_entity_key,only_visit_ghosts_once);
      }
      for ( std::vector<EntityProc>::const_iterator i = shared_change.begin() ; i != shared_change.end() ; ++i) {
          impl::VisitAuraClosureGeneral(*this,i->first,store_entity_key,only_visit_ghosts_once);
      }
      for ( std::set<EntityProc,EntityLess>::const_iterator i = send_closure.begin() ; i != send_closure.end() ; ++i) {
          impl::VisitAuraClosureGeneral(*this,i->first,store_entity_key,only_visit_ghosts_once);
      }

    std::set<EntityKey> & modified_ghosts = store_entity_key.entity_key_set;

    // The ghosted change list will become invalid
    ghosted_change.clear();

    std::vector<EntityProc> empty_add ;
    std::vector<EntityKey>  remove_modified_ghosts( modified_ghosts.begin() ,
                                                    modified_ghosts.end() );

    // Skip 'm_ghosting[0]' which is the shared subset.
    for ( std::vector<Ghosting*>::iterator
          ig = m_ghosting.begin() + 1; ig != m_ghosting.end() ; ++ig ) {
      // parallel synchronous:
      internal_change_ghosting( **ig , empty_add , remove_modified_ghosts );
    }
  }



  //------------------------------
  // Consistently change the owner on all processes.
  // 1) The local_change list is giving away ownership.
  // 2) The shared_change may or may not be receiving ownership

  {
    PartVector owned;
    owned.push_back(& meta.locally_owned_part());

    for ( std::vector<EntityProc>::iterator
          i = local_change.begin() ; i != local_change.end() ; ++i ) {
      // Giving ownership, change the parts first and then
      // the owner rank to pass the ownership test.
      Entity entity = i->first;

      internal_verify_and_change_entity_parts( entity , PartVector() , owned );

      const bool changed = this->set_parallel_owner_rank_but_not_comm_lists( entity, i->second );
      if (changed) {
        internal_change_owner_in_comm_data(entity_key(entity), i->second);
      }
    }

    for ( std::vector<EntityProc>::iterator
          i = shared_change.begin() ; i != shared_change.end() ; ++i ) {
      Entity entity = i->first;
      const bool changed = this->set_parallel_owner_rank_but_not_comm_lists( entity, i->second );
      if (changed) {
        internal_change_owner_in_comm_data(entity_key(entity), i->second);
      }
      if ( p_rank == i->second ) { // I receive ownership
          internal_verify_and_change_entity_parts( entity , owned , PartVector() );
      }
    }
  }


  //------------------------------
  // Send entities, along with their closure, to the new owner processes
  {
    std::ostringstream error_msg ;
    int error_count = 0 ;

    CommAll comm( p_comm );

    EntityVector unique_list_of_send_closure;
    unique_list_of_send_closure.reserve(send_closure.size());

    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      CommBuffer & buffer = comm.send_buffer( i->second );
      Entity entity = i->first;
      pack_entity_info(*this, buffer , entity );
      pack_field_values(*this, buffer , entity );

      if (unique_list_of_send_closure.empty() || entity_key(unique_list_of_send_closure.back()) != entity_key(entity)) {
        unique_list_of_send_closure.push_back(entity);
      }
    }

    BABBLE_STK_PARALLEL_COMM(p_comm, "          change_entity_owner calling allocate_buffers");
    comm.allocate_buffers( p_size / 4 );

    for ( std::set<EntityProc,EntityLess>::iterator
          i = send_closure.begin() ; i != send_closure.end() ; ++i ) {
      CommBuffer & buffer = comm.send_buffer( i->second );
      Entity entity = i->first;
      pack_entity_info(*this, buffer , entity );
      pack_field_values(*this, buffer , entity );
    }

    BABBLE_STK_PARALLEL_COMM(p_comm, "          change_entity_owner calling communicate");
    comm.communicate();

    for ( int p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = comm.recv_buffer(p);
      while ( buf.remaining() ) {
        PartVector parts ;
        std::vector<Relation> relations ;
        EntityKey key ;
        int owner = ~0u ;

        unpack_entity_info( buf, *this, key, owner, parts, relations );

        // Received entity information will be correct,
        // modulo the owned and shared parts

        remove( parts , meta.globally_shared_part() );

        if ( owner == p_rank ) {
          // Must have the locally_owned_part
          insert( parts , meta.locally_owned_part() );
        }
        else {
          // Must not have the locally_owned_part
          remove( parts , meta.locally_owned_part() );
        }

        std::pair<Entity ,bool> result = internal_create_entity( key );

        Entity entity = result.first;

        // The entity was copied and not created.

        internal_change_entity_parts( entity , parts , PartVector() );

        log_created_parallel_copy( entity );

        const bool changed = this->set_parallel_owner_rank_but_not_comm_lists( entity, owner );
        if (changed) {
          internal_change_owner_in_comm_data(entity_key(entity), owner);
        }

        internal_declare_relation( entity , relations );

        if ( ! unpack_field_values(*this, buf , entity , error_msg ) ) {
          ++error_count ;
        }
      }
    }

    all_reduce( p_comm , ReduceSum<1>( & error_count ) );
    ThrowErrorMsgIf( error_count, error_msg.str() );

    // Any entity that I sent and is not in an owned closure is deleted.
    // The owned closure will be effected by received entities, so can
    // only clean up after the newly owned entities have been received.
    // Destroy backwards so as not to invalidate closures in the process.

    {
      for ( EntityVector::reverse_iterator i = unique_list_of_send_closure.rbegin() ;
            i != unique_list_of_send_closure.rend() ;
            ++i) {
        if ( ! stk::mesh::impl::member_of_owned_closure(*this, *i , p_rank ) ) {
          ThrowRequireMsg( internal_destroy_entity( *i ),
                           "Failed to destroy entity " << identifier(*i) );
        }
      }
    }

    send_closure.clear(); // Has been invalidated
  }

#ifdef CHANGE_ENTITY_OWNER_UPDATES_COMM_INFO
  internal_apply_node_sharing(node_to_dependent_processors_map);
  internal_print_comm_map("A-APNS");
  internal_print_comm_list("A-APNS", *this);
#endif


}

//----------------------------------------------------------------------


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//Ghosting
//----------------------------------------------------------------------

//----------------------------------------------------------------------

Ghosting & BulkData::create_ghosting( const std::string & name )
{
    INCREMENT_MODIFICATION_COUNTER(PUBLIC, CREATE_GHOSTING);
    return internal_create_ghosting(name);
}

Ghosting & BulkData::internal_create_ghosting( const std::string & name )
{
    INCREMENT_MODIFICATION_COUNTER(INTERNAL, CREATE_GHOSTING);

  require_ok_to_modify();

  // Verify name is the same on all processors,
  // if not then throw an exception on all processors.
  if (parallel_size() > 1) {
    CommBroadcast bc( parallel() , 0 );

    if ( bc.parallel_rank() == 0 ) {
      bc.send_buffer().skip<char>( name.size() + 1 );
    }

    bc.allocate_buffer();

    if ( bc.parallel_rank() == 0 ) {
      bc.send_buffer().pack<char>( name.c_str() , name.size() + 1 );
    }

    bc.communicate();

    const char * const bc_name =
      reinterpret_cast<const char *>( bc.recv_buffer().buffer() );

    int error = 0 != std::strcmp( bc_name , name.c_str() );

    all_reduce( parallel() , ReduceMax<1>( & error ) );

    ThrowErrorMsgIf( error, "Parallel name inconsistency");
  }

  Ghosting * const g =
    new Ghosting( *this , name , m_ghosting.size() , m_sync_count );

  m_ghosting.push_back( g );

  if (m_ghost_parts.size() == 0) {
    ThrowRequireMsg(equal_case(std::string("shared"), name), "Expect shared to be the first ghosting created.");
    m_ghost_parts.push_back(&mesh_meta_data().globally_shared_part());
  }
  else if (m_ghost_parts.size() == 1) {
    ThrowRequireMsg(equal_case(std::string("shared_aura"), name), "Expect aura to be the second ghosting created.");
    Part & aura_part = mesh_meta_data().aura_part();
    aura_part.entity_membership_is_parallel_consistent(false);
    m_ghost_parts.push_back(&aura_part);
  }
  else {
    std::ostringstream oss;
    oss << "custom_ghosting_" << m_ghost_parts.size();
    std::string ghostPartName = stk::mesh::impl::convert_to_internal_name(oss.str());
    Part& ghost_part = mesh_meta_data().declare_part(ghostPartName);
    ghost_part.entity_membership_is_parallel_consistent(false);
    m_ghost_parts.push_back(&ghost_part);
  }

  ThrowRequireMsg(m_ghost_parts.size() == m_ghosting.size(), "m_ghost_parts.size()="<<m_ghost_parts.size()<<", must be same as m_ghosting.size()="<<m_ghosting.size());

  return *g ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void comm_recv_to_send(
  BulkData & mesh ,
  const std::set< EntityKey > & entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs ,
        std::set< EntityProc , EntityLess > & entitiesToGhostOntoOtherProcessors );

void comm_sync_send_recv(
  BulkData & mesh ,
  std::set< EntityProc , EntityLess > & new_send ,
  std::set< EntityKey > & new_recv );

} // namespace <>

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::destroy_ghosting( Ghosting& ghost_layer )
{
  INCREMENT_MODIFICATION_COUNTER(PUBLIC, DESTROY_GHOSTING);
  std::vector<EntityKey> receive_list;
  ghost_layer.receive_list(receive_list);
  internal_verify_inputs_and_change_ghosting(ghost_layer, std::vector<stk::mesh::EntityProc>(), receive_list);
}

//----------------------------------------------------------------------

void BulkData::destroy_all_ghosting()
{
  INCREMENT_MODIFICATION_COUNTER(PUBLIC, DESTROY_ALL_GHOSTING);

  Trace_("stk::mesh::BulkData::destroy_all_ghosting");

  require_ok_to_modify();

  // Clear Ghosting data

  for ( std::vector<Ghosting*>::iterator
        ig = m_ghosting.begin() ; ig != m_ghosting.end() ; ++ig ) {
    Ghosting & gh = **ig ;
    gh.m_sync_count = m_sync_count ;
  }

  // Iterate backwards so as not to invalidate a closure.

  for ( EntityCommListInfoVector::reverse_iterator
        i =  m_entity_comm_list.rbegin() ;
        i != m_entity_comm_list.rend() ; ++i) {

    if ( in_receive_ghost( i->key ) ) {
        entity_comm_map_clear( i->key );
        internal_destroy_entity( i->entity );
      i->key = EntityKey();
      i->entity_comm = NULL;
    }
    else {
        entity_comm_map_clear_ghosting(i->key);
      if ( entity_comm_map(i->key).empty() ) {
        i->key = EntityKey();
        i->entity_comm = NULL;
      }
    }
  }

  EntityCommListInfoVector::iterator i =
    std::remove_if( m_entity_comm_list.begin() ,
                    m_entity_comm_list.end() , IsInvalid() );

  m_entity_comm_list.erase( i , m_entity_comm_list.end() );
}

//----------------------------------------------------------------------

void BulkData::change_ghosting(
  Ghosting & ghosts ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<EntityKey> & remove_receive )
{
    INCREMENT_MODIFICATION_COUNTER(PUBLIC, CHANGE_GHOSTING);
    internal_verify_inputs_and_change_ghosting(ghosts, add_send, remove_receive);
}

void BulkData::internal_verify_inputs_and_change_ghosting(
  Ghosting & ghosts ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<EntityKey> & remove_receive )
{
  Trace_("stk::mesh::BulkData::change_ghosting");

  //----------------------------------------
  // Verify inputs:

  require_ok_to_modify();

  const bool ok_mesh  = & BulkData::get(ghosts) == this ;
  const bool ok_ghost = 1 < ghosts.ordinal();
  bool ok_add    = true ;
  bool ok_remove = true ;

  // Verify all 'add' are locally owned.

  for ( std::vector<EntityProc>::const_iterator
        i = add_send.begin() ; ok_add && i != add_send.end() ; ++i ) {
    ok_add = parallel_owner_rank(i->first) == parallel_rank();
  }

  // Verify all 'remove' are members of the ghosting.

  for ( std::vector<EntityKey>::const_iterator
        i = remove_receive.begin() ;
        ok_remove && i != remove_receive.end() ; ++i ) {
    ok_remove = in_receive_ghost( ghosts , *i );
  }

  int ok = ok_mesh && ok_ghost && ok_add && ok_remove ;

  all_reduce( parallel() , ReduceMin<1>( & ok ) );

  if ( 0 == ok ) {
    std::ostringstream msg ;
    msg << "For ghosts " << ghosts.name() << ", " ;
    if ( ! ok_mesh )  { msg << " : Mesh does not own this ghosting" ; }
    if ( ! ok_ghost ) { msg << " : Cannot modify this ghosting" ; }
    if ( ! ok_add ) {
      msg << " : Not owned add {" ;
      for ( std::vector<EntityProc>::const_iterator
            i = add_send.begin() ; i != add_send.end() ; ++i ) {
        if ( parallel_owner_rank(i->first) != parallel_rank() ) {
          msg << " " << identifier(i->first);
        }
      }
      msg << " }" ;
    }
    if ( ! ok_remove ) {
      msg << " : Not in ghost receive {" ;
      for ( std::vector<EntityKey>::const_iterator
            i = remove_receive.begin() ; i != remove_receive.end() ; ++i ) {
        if ( ! in_receive_ghost( ghosts , *i ) ) {
          msg << " " << i->id();
        }
      }
    }

    ThrowErrorMsg( msg.str() );
  }
  //----------------------------------------
  // Change the ghosting:

  internal_change_ghosting( ghosts , add_send , remove_receive );
}

//----------------------------------------------------------------------

void BulkData::ghost_entities_and_fields(Ghosting & ghosting, const std::set<EntityProc , EntityLess>& entitiesToGhostOntoOtherProcessors)
{
    const size_t record_entity_comm_size_before_changing_it = m_entity_comm_list.size();
    const int p_size = parallel_size() ;

    CommAll comm( m_parallel_machine );

    for ( int phase = 0; phase < 2; ++phase ) {
      for ( std::set< EntityProc , EntityLess >::iterator
              j = entitiesToGhostOntoOtherProcessors.begin(); j != entitiesToGhostOntoOtherProcessors.end() ; ++j ) {

        Entity entity = j->first;
        const int proc = j->second;

        if ( ! in_ghost( ghosting , entity_key(entity) , proc ) ) {
          // Not already being sent , must send it.
          CommBuffer & buf = comm.send_buffer( proc );
          buf.pack<unsigned>( entity_rank(entity) );
          pack_entity_info(*this, buf , entity );
          pack_field_values(*this, buf , entity );

          if (phase == 1) {
            entity_comm_map_insert(entity, EntityCommInfo(ghosting.ordinal(), proc));
            const EntityComm* entity_comm = m_entity_comm_map.entity_comm(entity_key(entity));
            EntityCommListInfo comm_info = {entity_key(entity), entity, parallel_owner_rank(entity), entity_comm};
            m_entity_comm_list.push_back( comm_info );
          }
        }
      }

      if (phase == 0) {
        comm.allocate_buffers( p_size / 4 );
      }
      else {
        comm.communicate();
      }
    }

    std::ostringstream error_msg ;
    int error_count = 0 ;
    OrdinalVector ordinal_scratch;
    PartVector part_scratch;

    const MetaData & meta = m_mesh_meta_data ;
    const unsigned rank_count = meta.entity_rank_count();

    for ( unsigned rank = 0 ; rank < rank_count ; ++rank ) {
      for ( int p = 0 ; p < p_size ; ++p ) {
        CommBuffer & buf = comm.recv_buffer(p);
        while ( buf.remaining() ) {
          // Only unpack if of the current entity rank.
          // If not the current entity rank, break the iteration
          // until a subsequent entity rank iteration.
          {
            unsigned this_rank = ~0u ;
            buf.peek<unsigned>( this_rank );

            if ( this_rank != rank ) break ;

            buf.unpack<unsigned>( this_rank );
          }

          PartVector parts ;
          std::vector<Relation> relations ;
          EntityKey key ;
          int owner = ~0u ;

          unpack_entity_info( buf, *this, key, owner, parts, relations );

          // Must not have the locally_owned_part or globally_shared_part

          remove( parts , meta.locally_owned_part() );
          remove( parts , meta.globally_shared_part() );

          if (owner != this->parallel_rank()) {
            // We will also add the entity to the part corresponding to the 'ghosts' ghosting.
            stk::mesh::Part& ghost_part = *m_ghost_parts[ghosting.ordinal()];
            insert( parts, ghost_part );
          }

          GhostReuseMap::iterator f_itr = m_ghost_reuse_map.find(key);
          const size_t use_this_offset = f_itr == m_ghost_reuse_map.end() ? 0 : f_itr->second;
          if (use_this_offset != 0) {
            m_ghost_reuse_map.erase(f_itr);
          }

          std::pair<Entity ,bool> result = internal_create_entity( key, use_this_offset );

          Entity entity = result.first;
          const bool created   = result.second ;

          require_entity_owner( entity , owner );

          internal_change_entity_parts( entity , parts , PartVector() );

          if ( created ) {
            log_created_parallel_copy( entity );
            this->set_parallel_owner_rank_but_not_comm_lists( entity, owner);
          }

          internal_declare_relation( entity , relations, ordinal_scratch, part_scratch );

          if ( ! unpack_field_values(*this, buf , entity , error_msg ) ) {
            ++error_count ;
          }

          const EntityCommInfo tmp( ghosting.ordinal() , owner );

          if ( entity_comm_map_insert(entity, tmp) ) {
            const EntityComm* entity_comm = m_entity_comm_map.entity_comm(entity_key(entity));
            EntityCommListInfo comm_info = {entity_key(entity), entity, parallel_owner_rank(entity), entity_comm};
            m_entity_comm_list.push_back( comm_info );
          }
        }
      }
    }

    if (parallel_size() > 1) {
      all_reduce( m_parallel_machine , ReduceSum<1>( & error_count ) );
    }

    ThrowErrorMsgIf( error_count, error_msg.str() );

    if ( record_entity_comm_size_before_changing_it < m_entity_comm_list.size() ) {
      // Added new ghosting entities to the list,
      // must now sort and merge.

      EntityCommListInfoVector::iterator i = m_entity_comm_list.begin();
      i += record_entity_comm_size_before_changing_it ;
      std::sort( i , m_entity_comm_list.end() );
      std::inplace_merge( m_entity_comm_list.begin() , i ,
                          m_entity_comm_list.end() );
      m_entity_comm_list.erase( std::unique( m_entity_comm_list.begin() , m_entity_comm_list.end() ) ,
                                m_entity_comm_list.end() );

      internal_sync_comm_list_owners();
    }
    ghosting.m_sync_count = m_sync_count ;
}

namespace impl {

struct StoreInEntityProcSet {
    StoreInEntityProcSet(
            BulkData & mesh_in,
            std::set<stk::mesh::EntityProc, stk::mesh::EntityLess> & set_in)
    :mesh(mesh_in)
    ,myset(set_in) { }

    void operator()(Entity entity) {
      myset.insert(stk::mesh::EntityProc(entity,proc));
    }

    BulkData & mesh;
    std::set<stk::mesh::EntityProc , stk::mesh::EntityLess> & myset;
    int proc;
};

struct OnlyGhosts  {
    OnlyGhosts(BulkData & mesh_in) : mesh(mesh_in) {}
    bool operator()(Entity entity) {
        const bool isValid = mesh.is_valid(entity);
        const bool iDoNotOwnEntity = proc != mesh.parallel_owner_rank(entity);
        const bool entityIsShared = mesh.in_shared( mesh.entity_key(entity) , proc );
        return (isValid && iDoNotOwnEntity && !entityIsShared);
    }
    BulkData & mesh;
    int proc;
};
} // namespace impl

void BulkData::internal_change_ghosting(
  Ghosting & ghosting ,
  const std::vector<EntityProc> & add_send ,
  const std::vector<EntityKey> & remove_receive,
  bool is_full_regen)
{
  Trace_("stk::mesh::BulkData::internal_change_ghosting");

  BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "      entered internal_change_ghosting");

  INCREMENT_MODIFICATION_COUNTER(INTERNAL, CHANGE_GHOSTING);

  //------------------------------------
  // Copy ghosting lists into more efficiently edited container.
  // The send and receive lists must be in entity rank-order.

  std::set<EntityProc , EntityLess> entitiesToGhostOntoOtherProcessors(EntityLess(*this));
  std::set<EntityKey>               entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs;

  //------------------------------------
  // Insert the current ghost receives and then remove from that list.

  // This if-check is an optimization; if doing a full regen
  // then we are removing all ghosting information and new_recv should
  // be left empty.
  if ( !is_full_regen ) {

    // Iterate over all entities with communication information, adding
    // the entity if it's a ghost on this process. new_recv will contain
    // all ghosts on this process by the end of the loop.
    for ( EntityCommListInfoVector::const_iterator
          i = comm_list().begin() ; i != comm_list().end() ; ++i ) {
      if ( in_receive_ghost( ghosting , i->key ) ) {
        entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.insert( i->key );
      }
    }

    // Remove any entities that are in the remove list.

    for ( std::vector<EntityKey>::const_iterator
          i = remove_receive.begin() ; i != remove_receive.end() ; ++i ) {
      entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.erase( *i );
    }

    // Keep the closure of the remaining received ghosts.
    // Working from highest-to-lowest key (rank entity type)
    // results in insertion of the closure because
    // inserted entities will get looped over after they are inserted.

    // Insertion will not invalidate the associative container's iterator.

    for ( std::set<EntityKey>::reverse_iterator
          i = entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.rbegin() ; i != entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.rend() ; ++i) {
      const unsigned erank = i->rank();

      Entity e = get_entity(*i); // Could be performance issue? Not if you're just doing full regens

      for (EntityRank irank = stk::topology::BEGIN_RANK; irank < erank; ++irank) {
        Entity const *rels_i = begin(e, irank);
        Entity const *rels_e = end(e, irank);
        for (; rels_i != rels_e; ++rels_i)
        {
          if ( is_valid(*rels_i) && in_receive_ghost( ghosting , entity_key(*rels_i) ) )
          {
            entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.insert( entity_key(*rels_i) );
          }
        }
      }
    }
  }

  //  Initialize the new_send from the new_recv
  comm_recv_to_send( *this , entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs , entitiesToGhostOntoOtherProcessors );

  //------------------------------------
  // Add the specified entities and their closure to the send ghosting

  impl::StoreInEntityProcSet sieps(*this,entitiesToGhostOntoOtherProcessors);
  impl::OnlyGhosts og(*this);
  for ( std::vector< EntityProc >::const_iterator
        i = add_send.begin() ; i != add_send.end() ; ++i ) {
      og.proc = i->second;
      sieps.proc = i->second;
      impl::VisitClosureGeneral(*this,i->first,sieps,og);
  }

  // Synchronize the send and receive list.
  // If the send list contains a not-owned entity
  // inform the owner and receiver to add that entity
  // to their ghost send and receive lists.

  comm_sync_send_recv( *this , entitiesToGhostOntoOtherProcessors , entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs );

  // The new_send list is now parallel complete and parallel accurate
  // The new_recv has those ghost entities that are to be kept.
  //------------------------------------
  // Remove the ghost entities that will not remain.
  // If the last reference to the receive ghost entity then delete it.

  bool removed = false ;

  for ( EntityCommListInfoVector::reverse_iterator
        i = m_entity_comm_list.rbegin() ; i != m_entity_comm_list.rend() ; ++i) {

    const bool is_owner = i->owner == m_parallel_rank ;
    const bool remove_recv = ( ! is_owner ) &&
                             0 == entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.count(i->key);

    if ( is_owner ) {
      // Is owner, potentially removing ghost-sends
      // Have to make a copy

      std::vector<EntityCommInfo> comm_ghost ;
      const PairIterEntityComm ec = entity_comm_map(i->key, ghosting);
      comm_ghost.assign( ec.first , ec.second );

      for ( ; ! comm_ghost.empty() ; comm_ghost.pop_back() ) {
        const EntityCommInfo tmp = comm_ghost.back();

        if ( 0 == entitiesToGhostOntoOtherProcessors.count( EntityProc( i->entity , tmp.proc ) ) ) {
          entity_comm_map_erase(i->key, tmp);
        }
      }
    }
    else if ( remove_recv ) {
        entity_comm_map_erase(i->key, ghosting);
    }

    if ( entity_comm_map(i->key).empty() ) {
      removed = true ;
      i->key = EntityKey(); // No longer communicated
      if ( remove_recv ) {
        ThrowRequireMsg( internal_destroy_entity( i->entity, remove_recv ),
                         " FAILED attempt to destroy entity: " << identifier(i->entity) );
      }
    }
  }

  if ( removed ) {
    EntityCommListInfoVector::iterator i =
      std::remove_if( m_entity_comm_list.begin() ,
                      m_entity_comm_list.end() , IsInvalid() );
    m_entity_comm_list.erase( i , m_entity_comm_list.end() );
  }

  //------------------------------------
  // Push newly ghosted entities to the receivers and update the comm list.
  // Unpacking must proceed in entity-rank order so that higher ranking
  // entities that have relations to lower ranking entities will have
  // the lower ranking entities unpacked first.  The higher and lower
  // ranking entities may be owned by different processes,
  // as such unpacking must be performed in rank order.

  ghost_entities_and_fields(ghosting, entitiesToGhostOntoOtherProcessors);
}

//----------------------------------------------------------------------


namespace {

// Fill a new send list from the receive list.
void comm_recv_to_send(
  BulkData & mesh ,
  const std::set< EntityKey > & entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs,
        std::set< EntityProc , EntityLess > & entitiesToGhostOntoOtherProcessors )
{
  const int parallel_size = mesh.parallel_size();

  CommAll all( mesh.parallel() );

  // For all entity keys in new recv, send the entity key to the owning proc
  for ( int phase = 0; phase < 2; ++phase) {
    for ( std::set<EntityKey>::const_iterator
            i = entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.begin() ; i != entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.end() ; ++i ) {
      Entity e = mesh.get_entity(*i); // Performance issue? Not if you're just doing full regens
      const int owner = mesh.parallel_owner_rank(e);
      const EntityKey key = mesh.entity_key(e);
      all.send_buffer( owner ).pack<EntityKey>( key );
    }
    if (phase == 0) { //allocation phase
      BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "          comm_recv_to_send calling allocate_buffers");
      all.allocate_buffers( parallel_size / 4 , false /* Not symmetric */ );
    }
    else { //communication phase
      BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "          comm_recv_to_send calling communicate");
      all.communicate();
    }
  }

  // Insert into entitiesToGhostOntoOtherProcessors that there is an entity on another proc
  for ( int proc_rank = 0 ; proc_rank < parallel_size ; ++proc_rank ) {
    CommBuffer & buf = all.recv_buffer(proc_rank);
    while ( buf.remaining() ) {
      EntityKey key ;
      buf.unpack<EntityKey>( key );
      EntityProc tmp( mesh.get_entity( key ) , proc_rank );
      entitiesToGhostOntoOtherProcessors.insert( tmp );
    }
  }
}

// Synchronize the send list to the receive list.

void comm_sync_send_recv(
  BulkData & mesh ,
  std::set< EntityProc , EntityLess > & entitiesToGhostOntoOtherProcessors ,
  std::set< EntityKey > & entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs )
{
  const int parallel_rank = mesh.parallel_rank();
  const int parallel_size = mesh.parallel_size();

  CommAll all( mesh.parallel() );

  // Communication sizing:

  for ( std::set< EntityProc , EntityLess >::iterator
        i = entitiesToGhostOntoOtherProcessors.begin() ; i != entitiesToGhostOntoOtherProcessors.end() ; ++i ) {
    const int owner = mesh.parallel_owner_rank(i->first);
    all.send_buffer( i->second ).skip<EntityKey>(1).skip<int>(1);
    if ( owner != parallel_rank ) {
      all.send_buffer( owner ).skip<EntityKey>(1).skip<int>(1);
    }
  }

  BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "          comm_sync_send_recv calling allocate_buffers");
  all.allocate_buffers( parallel_size / 4 , false /* Not symmetric */ );

  // Loop thru all entities in entitiesToGhostOntoOtherProcessors, send the entity key to the sharing/ghosting proc
  // Also, if the owner of the entity is NOT me, also send the entity key to the owing proc

  // Communication packing (with message content comments):
  for ( std::set< EntityProc , EntityLess >::iterator
        i = entitiesToGhostOntoOtherProcessors.begin() ; i != entitiesToGhostOntoOtherProcessors.end() ; ) {
    const int owner = mesh.parallel_owner_rank(i->first);

    // Inform receiver of ghosting, the receiver does not own
    // and does not share this entity.
    // The ghost either already exists or is a to-be-done new ghost.
    // This status will be resolved on the final communication pass
    // when new ghosts are packed and sent.

    const EntityKey entity_key = mesh.entity_key(i->first);
    const int proc = i->second;

    all.send_buffer( i->second ).pack(entity_key).pack(proc);

    if ( owner != parallel_rank ) {
      // I am not the owner of this entity.
      // Inform the owner of this ghosting need.
      all.send_buffer( owner ).pack(entity_key).pack(proc);

      // Erase it from my processor's ghosting responsibility:
      // The iterator passed to the erase method will be invalidated.
      std::set< EntityProc , EntityLess >::iterator jrem = i ; ++i ;
      entitiesToGhostOntoOtherProcessors.erase( jrem );
    }
    else {
      ++i ;
    }
  }

  BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "          comm_sync_send_recv calling communicate");
  all.communicate();

  // Loop thru all the buffers, and insert ghosting request for entity e to other proc
  // if the proc sending me the data is me, then insert into new_recv.
  // Communication unpacking:
  for ( int p = 0 ; p < parallel_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer(p);
    while ( buf.remaining() ) {

      EntityKey entity_key;
      int proc = 0;

      buf.unpack(entity_key).unpack(proc);

      Entity const e = mesh.get_entity( entity_key );

      if ( parallel_rank != proc ) {
        //  Receiving a ghosting need for an entity I own.
        //  Add it to my send list.
        ThrowRequireMsg(mesh.is_valid(e),
            "Unknown entity key: " <<
            MetaData::get(mesh).entity_rank_name(entity_key.rank()) <<
            "[" << entity_key.id() << "]");
        EntityProc tmp( e , proc );
        entitiesToGhostOntoOtherProcessors.insert( tmp );
      }
      else if ( mesh.is_valid(e) ) {
        //  I am the receiver for this ghost.
        //  If I already have it add it to the receive list,
        //  otherwise don't worry about it - I will receive
        //  it in the final new-ghosting communication.
        entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs.insert( mesh.entity_key(e) );
      }
    }
  }
}

void insert_upward_relations(const BulkData& bulk_data, Entity rel_entity,
                             const EntityRank rank_of_orig_entity,
                             const int my_rank,
                             const int share_proc,
                             std::vector<EntityProc>& send)
{
  EntityRank rel_entity_rank = bulk_data.entity_rank(rel_entity);
  ThrowAssert(rel_entity_rank > rank_of_orig_entity);

  // If related entity is higher rank, I own it, and it is not
  // already shared by proc, ghost it to the sharing processor.
  if ( bulk_data.parallel_owner_rank(rel_entity) == my_rank &&
       ! bulk_data.in_shared( bulk_data.entity_key(rel_entity) , share_proc ) ) {

    EntityProc entry( rel_entity , share_proc );
    send.push_back( entry );

    // There may be even higher-ranking entities that need to be ghosted, so we must recurse
    const EntityRank end_rank = static_cast<EntityRank>(bulk_data.mesh_meta_data().entity_rank_count());
    EntityVector temp_entities;
    Entity const* rels = NULL;
    int num_rels = 0;
    for (EntityRank irank = static_cast<EntityRank>(rel_entity_rank + 1); irank < end_rank; ++irank)
    {
      if (bulk_data.connectivity_map().valid(rel_entity_rank, irank)) {
        num_rels = bulk_data.num_connectivity(rel_entity, irank);
        rels     = bulk_data.begin(rel_entity, irank);
      }
      else {
        num_rels = get_connectivity(bulk_data, rel_entity, irank, temp_entities);
        rels     = &*temp_entities.begin();
      }

      for (int r = 0; r < num_rels; ++r)
      {
        Entity const rel_of_rel_entity = rels[r];
        if (bulk_data.is_valid(rel_of_rel_entity)) {
          insert_upward_relations(bulk_data, rel_of_rel_entity, rel_entity_rank, my_rank, share_proc, send);
        }
      }
    }
  }
}

} // namespace <>

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::internal_regenerate_aura()
{
  Trace_("stk::mesh::BulkData::internal_regenerate_shared_aura");
  BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "  entered internal_regenerate_shared_aura");

  require_ok_to_modify();

  std::vector<EntityProc> send ;

  const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());

  // Iterate over all entities with communication info, get the sharing
  // comm info for each entity, and ensure that upwardly related
  // entities to the shared entity are ghosted on the sharing proc.
  EntityVector temp_entities;
  Entity const* rels = NULL;
  int num_rels = 0;
  for ( EntityCommListInfoVector::const_iterator
      i = comm_list().begin() ; i != comm_list().end() ; ++i ) {

    const EntityRank erank = static_cast<EntityRank>(i->key.rank());

    const PairIterEntityComm aura = entity_comm_map_shared(i->key);

    for ( size_t j = 0 ; j < aura.size() ; ++j ) {

      const int share_proc = aura[j].proc ;

      for (EntityRank k_rank = static_cast<EntityRank>(erank + 1); k_rank < end_rank; ++k_rank)
      {
        if (connectivity_map().valid(erank, k_rank)) {
          num_rels = num_connectivity(i->entity, k_rank);
          rels     = begin(i->entity, k_rank);
        }
        else {
          num_rels = get_connectivity(*this, i->entity, k_rank, temp_entities);
          rels     = &*temp_entities.begin();
        }

        for (int r = 0; r < num_rels; ++r)
        {
          if (is_valid(rels[r])) {
            insert_upward_relations(*this, rels[r], erank, m_parallel_rank, share_proc, send);
          }
        }
      }
    }
  }

  // Add new aura, remove all of the old aura.
  // The change_ghosting figures out what to actually delete and add.
  internal_change_ghosting( aura_ghosting() , send , std::vector<EntityKey>(), true /*full regen*/ );

  // BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "  exiting internal_regenerate_shared_aura");

}


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//EndSync
//----------------------------------------------------------------------


//----------------------------------------------------------------------

int BulkData::determine_new_owner( Entity entity ) const
{
  // We will decide the new owner by looking at all the processes sharing
  // this entity. The new owner will be the sharing process with lowest rank.

  // The local process is a candidate only if the entity is not destroyed.
  int new_owner = is_valid(entity) ? m_parallel_rank : ~0u;

  for ( PairIterEntityComm
        share = entity_comm_map_shared(entity_key(entity)); ! share.empty() ; ++share ) {
    if ( share->proc < m_parallel_size &&
         ( share->proc < new_owner || m_parallel_size <= new_owner ) ) {
      new_owner = share->proc ;
    }
  }

  return new_owner ;
}

//----------------------------------------------------------------------

namespace {

// A method for quickly finding an entity within a comm list
EntityCommListInfo find_entity(const BulkData& mesh,
                               const EntityCommListInfoVector& entities,
                               const EntityKey& key)
{
  EntityCommListInfoVector::const_iterator lb_itr = std::lower_bound(entities.begin(), entities.end(), key);
  ThrowAssertMsg(lb_itr != entities.end() && lb_itr->key == key,
                 "Cannot find id: " << key.id() << " in comm-list" );
  return *lb_itr;
}

struct EntityParallelState {
  int                 from_proc;
  EntityState         state;
  EntityCommListInfo  comm_info;
  const BulkData* mesh;

  bool operator<(const EntityParallelState& rhs) const
  { return EntityLess(*mesh)(comm_info.entity, rhs.comm_info.entity); }
};


bool pack_entity_modification( const BulkData & mesh ,
                               const bool packShared ,
                               CommAll & comm )
{
  bool flag = false;
  bool packGhosted = packShared == false;

  const EntityCommListInfoVector & entityCommList = mesh.comm_list();

  for ( EntityCommListInfoVector::const_iterator
        i = entityCommList.begin() ; i != entityCommList.end() ; ++i ) {

    Entity entity = i->entity;
    EntityState status = mesh.is_valid(entity) ? mesh.state(entity) : Deleted;

    if ( status == Modified || status == Deleted ) {

      for ( PairIterEntityComm ec = mesh.entity_comm_map(i->key); ! ec.empty() ; ++ec )
      {
        if ( ( packGhosted && ec->ghost_id > 0 ) || ( packShared && ec->ghost_id == 0 ) )
        {
          comm.send_buffer( ec->proc )
              .pack<EntityKey>( i->key )
              .pack<EntityState>( status );

          flag = true ;
        }
      }
    }
  }

  return flag ;
}

void communicate_entity_modification( const BulkData & mesh ,
                                      const bool shared ,
                                      std::vector<EntityParallelState > & data )
{
  BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "      entered communicate_entity_modification");
  CommAll comm( mesh.parallel() );
  const int p_size = comm.parallel_size();

  // Sizing send buffers:
  const bool local_mod = pack_entity_modification( mesh , shared , comm );

  // Allocation of send and receive buffers:
  BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "          communicate_entity_modification calling allocate_buffers")
  const bool global_mod = comm.allocate_buffers( comm.parallel_size() / 4 , false , local_mod );

  if ( global_mod )
  {
    const EntityCommListInfoVector & entityCommList = mesh.comm_list();

    // Packing send buffers:
    pack_entity_modification( mesh , shared , comm );

    BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "          communicate_entity_modification calling communicate")
    comm.communicate();

    for ( int procNumber = 0 ; procNumber < p_size ; ++procNumber ) {
      CommBuffer & buf = comm.recv_buffer( procNumber );
      EntityKey key;
      EntityState state;

      while ( buf.remaining() ) {

        buf.unpack<EntityKey>( key )
           .unpack<EntityState>( state );

        // search through entity_comm, should only receive info on entities
        // that are communicated.
        EntityCommListInfo info = find_entity(mesh, entityCommList, key);
        EntityParallelState parallel_state = {procNumber, state, info, &mesh};
        data.push_back( parallel_state );
      }
    }
  }

  std::sort( data.begin() , data.end() );

  // BABBLE_STK_PARALLEL_COMM(mesh.parallel(), "     exiting communicate_entity_modification");
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

// Postconditions:
//  * DistributedIndex is updated based on entity creation/deletions in the
//    last modification cycle.
//  * Comm lists for shared entities are up-to-date.
//  * shared_new contains all entities that were modified/created on a
//    different process

void BulkData::internal_update_distributed_index(
        stk::mesh::EntityRank entityRank, std::vector<Entity> & shared_new )
{
  Trace_("stk::mesh::BulkData::internal_update_distributed_index");
  std::vector<EntityKey> entity_keys;
  //resolve_edge_sharing(*this, entity_keys);
  resolve_entity_sharing(entityRank, entity_keys);

  shared_new.clear();
  shared_new.resize(entity_keys.size());
  for (size_t i=0;i<entity_keys.size();i++)
  {
      shared_new[i] = get_entity(entity_keys[i]);
  }
}

void BulkData::fillLocallyCreatedOrModifiedEntities(parallel::DistributedIndex::KeyTypeVector &local_created_or_modified)
{
  size_t num_created_or_modified = 0;

  for ( impl::EntityRepository::const_iterator
        i = m_entity_repo.begin() ; i != m_entity_repo.end() ; ++i )
  {
    Entity entity = i->second ;

    if ( state(entity) != Unchanged &&  is_entity_marked(entity)==BulkData::NOT_MARKED &&
         in_owned_closure( *this, entity , m_parallel_rank ) )
    {
      // Has been changed and is in owned closure, may be shared
      ++num_created_or_modified;
    }
  }

  local_created_or_modified.reserve(num_created_or_modified);

  for ( impl::EntityRepository::const_iterator
        i = m_entity_repo.begin() ; i != m_entity_repo.end() ; ++i )
  {
    Entity entity = i->second ;

    if ( state(entity) != Unchanged &&  is_entity_marked(entity)==BulkData::NOT_MARKED &&
         in_owned_closure(*this, entity , m_parallel_rank ) )
    {
      // Has been changed and is in owned closure, may be shared
      local_created_or_modified.push_back( entity_key(entity) );
    }
  }
}

void BulkData::internal_update_distributed_index(
        std::vector<Entity> & shared_new )
{
    Trace_("stk::mesh::BulkData::internal_update_distributed_index");
    BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "      entered internal_update_distributed_index");

    // Scoop up a list of all nodes that have had their sharing information
    // provided directly by the end-user.  These need special handling so
    // that DistributedIndex isn't used to query sharing information.
    //
    std::vector<EntityKey> shared_nodes;
    impl::gather_shared_nodes(*this, shared_nodes);

    std::vector<shared_entity_type> shared_edges;
    impl::markEntitiesForResolvingSharingInfoUsingNodes(*this, stk::topology::EDGE_RANK, shared_edges);

    std::vector<shared_entity_type> shared_faces;
    impl::markEntitiesForResolvingSharingInfoUsingNodes(*this, stk::topology::FACE_RANK, shared_faces);

    // Generate a list of all entities that have been created or modified
    // locally, who must have their sharing information in the comm lists
    // updated by making queries to DistributedIndex.  This list *does not*
    // include nodes that have had their sharing information explicitly
    // provided or edges and faces that will have their sharing information
    // resolved through node sharing.
    //
    parallel::DistributedIndex::KeyTypeVector local_created_or_modified;
    fillLocallyCreatedOrModifiedEntities(local_created_or_modified);

    // Update distributed index's ownership and sharing information with
    // locally-owned/shared entities that have been created or modified.
    // We need to manually insert the explicitly-shared nodes into this list
    // so that DistributedIndex is notified about IDs that are now in use
    // (since it manages available IDs).  These nodes will not have their
    // sharing info updated through queries to DistributedIndex.
    //
    parallel::DistributedIndex::KeyTypeVector keys_requiring_update(local_created_or_modified);
    keys_requiring_update.reserve(keys_requiring_update.size() + shared_nodes.size() + shared_edges.size() + shared_faces.size() );
    keys_requiring_update.insert(keys_requiring_update.end(), shared_nodes.begin(), shared_nodes.end());
    for (size_t i = 0; i < shared_edges.size(); ++i)
    {
      keys_requiring_update.push_back(shared_edges[i].global_key);
    }
    for (size_t i = 0; i < shared_faces.size(); ++i)
    {
      keys_requiring_update.push_back(shared_faces[i].global_key);
    }

    parallel::DistributedIndex::KeyTypeVector::const_iterator begin = keys_requiring_update.begin();
    parallel::DistributedIndex::KeyTypeVector::const_iterator end = keys_requiring_update.end();
    m_entities_index.update_keys( begin, end );

    // Update the sharing and ownership information in the comm maps for
    // locally-created or modified entities through queries to DistributedIndex.
    //

    std::vector<EntityKey> entity_keys;  // List of shared and modified entities
#define UPDATE_DIST_IDX_USE_QUERY_TO_USAGE
#ifdef UPDATE_DIST_IDX_USE_QUERY_TO_USAGE
    if (parallel_size() > 1)
    {
        // Retrieve data regarding which processes use the local_created_or_modified
        // including this process.
        parallel::DistributedIndex::KeyProcVector global_created_or_modified ;
        m_entities_index.query_to_usage( local_created_or_modified ,
                                 global_created_or_modified );

        //------------------------------
        // Take the usage data and update the sharing comm lists
        {
            Entity entity = Entity();

            // Iterate over all global modifications to this entity, this vector is
            // sorted, so we're guaranteed that all modifications to a particular
            // entities will be adjacent in this vector.

            for ( parallel::DistributedIndex::KeyProcVector::iterator
                    i =  global_created_or_modified.begin() ;
                    i != global_created_or_modified.end() ; ++i )
            {
                EntityKey key( static_cast<EntityKey::entity_key_t>(i->first) );
                int modifying_proc = i->second;

                if ( m_parallel_rank != modifying_proc )
                {
                    // Another process also created or updated this entity.
                    // Only want to look up entities at most once
                    if ( !is_valid(entity) || entity_key(entity) != key )
                    {
                        // Have not looked this entity up by key
                        entity = get_entity( key );
                        entity_keys.push_back(key);
                    }

                    // Add the other_process to the entity's sharing info.
                    entity_comm_map_insert(entity, EntityCommInfo( stk::mesh::BulkData::SHARED, modifying_proc ));
                }
            }
        }
    }
#endif

    update_shared_entities_global_ids( shared_edges );
    update_shared_entities_global_ids( shared_faces );

    // Add the explicitly-shared nodes to our list of shared and modified entities.
    //
    entity_keys.reserve(entity_keys.size() + shared_nodes.size());
    entity_keys.insert(entity_keys.end(), shared_nodes.begin(), shared_nodes.end());

    // Add the edges that had their sharing information updated through
    // node sharing to our list of shared and modified entities.
    //
    for (size_t i=0; i<shared_edges.size(); ++i)
    {
        Entity entity = get_entity(shared_edges[i].global_key);
        if ( is_entity_marked(entity) == BulkData::IS_SHARED )
        {
            entity_keys.push_back(shared_edges[i].global_key);
        }
    }

    for (size_t i=0; i<shared_faces.size(); ++i)
    {
        Entity entity = get_entity(shared_faces[i].global_key);
        if ( is_entity_marked(entity) == BulkData::IS_SHARED )
        {
            entity_keys.push_back(shared_faces[i].global_key);
        }
    }

    // Reset our marking array for all entities now that all sharing information
    // has been properly resolved.
    //
    std::fill(m_mark_entity.begin(), m_mark_entity.end(), static_cast<int>(BulkData::NOT_MARKED));

    shared_new.clear();
    shared_new.resize(entity_keys.size());
    for (size_t i=0; i<entity_keys.size(); ++i)
    {
        shared_new[i] = get_entity(entity_keys[i]);
    }
}



//----------------------------------------------------------------------

void BulkData::internal_establish_new_owner(stk::mesh::Entity entity)
{
    const int new_owner = determine_new_owner(entity);

    const bool changed = this->set_parallel_owner_rank_but_not_comm_lists(entity, new_owner);
    if(changed)
    {
        internal_change_owner_in_comm_data(entity_key(entity), new_owner);
    }
    set_synchronized_count(entity, m_sync_count);
}

void BulkData::internal_update_parts_for_shared_entity(stk::mesh::Entity entity, const bool is_entity_shared, const bool did_i_just_become_owner)
{
    PartVector parts_to_add_entity_to , parts_to_remove_entity_from ;

    if ( !is_entity_shared ) {
      parts_to_remove_entity_from.push_back(& m_mesh_meta_data.globally_shared_part());
    }

    if ( did_i_just_become_owner ) {
      parts_to_add_entity_to.push_back( & m_mesh_meta_data.locally_owned_part() );
    }

    if ( ! parts_to_add_entity_to.empty() || ! parts_to_remove_entity_from.empty() ) {
      internal_change_entity_parts( entity , parts_to_add_entity_to , parts_to_remove_entity_from );
    }
}

// Resolve modifications for shared entities:
// If not locally destroyed and remotely modified
// then set to locally modified.
// If remotely destroyed then determine the new owner.
//
// Post condition:
//  Shared entities are in-sync with respect to modification state.
//  Shared communication lists are updated to reflect all deletions.
//  Ownership has been re-assigned as necessary for deletion
//  of shared entities.

void BulkData::internal_resolve_shared_modify_delete()
{
  Trace_("stk::mesh::BulkData::internal_resolve_shared_modify_delete");
  BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "  entered internal_resolve_shared_modify_delete");

  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

  stk::mesh::impl::delete_shared_entities_which_are_no_longer_in_owned_closure( *this );

  std::vector< EntityParallelState > remotely_modified_shared_entities ;

  // Communicate entity modification state for shared entities
  // the resulting vector is sorted by entity and process.
  const bool communicate_shared = true ;
  communicate_entity_modification( *this , communicate_shared , remotely_modified_shared_entities );

  // We iterate backwards over remote_mod to ensure that we hit the
  // higher-ranking entities first.
  for ( std::vector<EntityParallelState>::reverse_iterator
        i = remotely_modified_shared_entities.rbegin(); i != remotely_modified_shared_entities.rend() ; ) {

    Entity entity                = i->comm_info.entity;
    EntityKey key                = i->comm_info.key;
    int owner                    = i->comm_info.owner;
    const bool locally_destroyed = !is_valid(entity);
    bool remote_owner_destroyed  = false;

    // Iterate over all of this entity's remote changes
    for ( ; i != remotely_modified_shared_entities.rend() && i->comm_info.entity == entity ; ++i ) {

      const int remote_proc    = i->from_proc ;
      const bool remotely_destroyed = Deleted == i->state ;

      // When a shared entity is remotely modified or destroyed
      // then the local copy is also modified.  This modification
      // status is applied to all related higher ranking entities.

      if ( ! locally_destroyed ) {
        this->mark_entity_and_upward_related_entities_as_modified( entity );
      }

      // A shared entity is being deleted on the remote process.
      // Remove it from the sharing communication list.
      // Ownership changes are processed later, but we'll need
      // to know if the remote owner destroyed the entity in order
      // to correctly resolve ownership (it is not sufficient to just
      // look at the comm list of the entity since there is no
      // guarantee that the comm list is correct or up-to-date).

      if ( remotely_destroyed ) {
        entity_comm_map_erase( key, EntityCommInfo(SHARED,remote_proc) );

        // check if owner is destroying
        if ( owner == remote_proc ) {
          remote_owner_destroyed = true ;
        }
      }
    }

    // Have now processed all remote changes knowledge for this entity.

    if(!locally_destroyed)
    {
        const bool am_i_old_local_owner = m_parallel_rank == parallel_owner_rank(entity);

        if ( remote_owner_destroyed ) {
            internal_establish_new_owner(entity);
        }

        const bool am_i_new_local_owner = m_parallel_rank == parallel_owner_rank(entity);
        const bool did_i_just_become_owner = ( ! am_i_old_local_owner && am_i_new_local_owner );

        const bool is_entity_shared = !entity_comm_map_shared(key).empty();
        internal_update_parts_for_shared_entity(entity, is_entity_shared, did_i_just_become_owner);
    }
  } // remote mod loop

  // Erase all sharing communication lists for Destroyed entities:
  for ( EntityCommListInfoVector::const_reverse_iterator
        i = comm_list().rbegin() ; i != comm_list().rend() ; ++i) {
    if ( !is_valid(i->entity) ) {
        entity_comm_map_erase( i->key, shared_ghosting() );
    }
  }
}


//----------------------------------------------------------------------
// Resolve modifications for ghosted entities:
// If a ghosted entity is modified or destroyed on the owning
// process then the ghosted entity must be destroyed.
//
// Post condition:
//  Ghosted entities of modified or deleted entities are destroyed.
//  Ghosted communication lists are cleared to reflect all deletions.

void BulkData::internal_resolve_ghosted_modify_delete()
{
  Trace_("stk::mesh::BulkData::internal_resolve_ghosted_modify_delete");
  BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "  entered internal_resolve_ghosted_modify_delete");


  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");
  // Resolve modifications for ghosted entities:

  std::vector<EntityParallelState > remotely_modified_ghosted_entities ;

  // Communicate entity modification state for ghost entities
  const bool communicate_shared = false ;
  communicate_entity_modification( *this , communicate_shared , remotely_modified_ghosted_entities );

  const size_t ghosting_count = m_ghosting.size();

  std::vector< int > ghosting_change_flags( ghosting_count , 0 );

  // We iterate backwards over remote_mod to ensure that we hit the
  // higher-ranking entities first. This is important because higher-ranking
  // entities like element must be deleted before the nodes they have are
  // deleted.
  for ( std::vector<EntityParallelState>::reverse_iterator
        i = remotely_modified_ghosted_entities.rbegin(); i != remotely_modified_ghosted_entities.rend() ; ++i ) {
    Entity entity                 = i->comm_info.entity;
    const EntityKey key           = i->comm_info.key;
    const int      remote_proc    = i->from_proc;
    const bool     local_owner    = i->comm_info.owner == m_parallel_rank ;
    const bool remotely_destroyed = Deleted == i->state ;
    const bool locally_destroyed  = !is_valid(entity);

    if ( local_owner ) { // Sending to 'remote_proc' for ghosting

      if ( remotely_destroyed ) {

        // remove from ghost-send list

         // j=2, j=1,
        for ( size_t j = ghosting_count ; --j ; ) {
          if ( entity_comm_map_erase( key, EntityCommInfo( j , remote_proc ) ) ) {
            ghosting_change_flags[ j ] = true ;
          }
        }
      }

      // Remotely modified ghosts are ignored

    }
    else { // Receiving from 'remote_proc' for ghosting

      // Owner modified or destroyed, must locally destroy.

        // Manoj: could be error here. I think the loop below
        // should loop thru all ghostings, m_ghosting. But effect might be same.
        // would only impact custom ghostings.
      for ( PairIterEntityComm ec = entity_comm_map(key); !ec.empty() ; ++ec ) {
        ghosting_change_flags[ ec->ghost_id ] = true ;
      }

      // This is a receive ghost so the only communication information
      // is the ghosting information, can clear it all out.
      //entity_comm_map_clear(key);
      entity_comm_map_erase(key, aura_ghosting());

      if ( ! locally_destroyed && ! in_owned_closure(*this, entity , m_parallel_rank ) )
      {
          // If mesh modification causes a ghost entity to become
          // a member of an owned-closure then do not automatically
          // destroy it.  The new sharing status will be resolved
          // in 'internal_resolve_parallel_create'.
          const bool destroy_entity_successful = internal_destroy_entity(entity, true);
          ThrowRequireMsg(destroy_entity_successful,
                          "Could not destroy ghost entity " << identifier(entity));
      }
    }
  } // end loop on remote mod

  // Erase all ghosting communication lists for:
  // 1) Destroyed entities.
  // 2) Owned and modified entities.

  for ( EntityCommListInfoVector::const_reverse_iterator
        i = comm_list().rbegin() ; i != comm_list().rend() ; ++i) {

    Entity entity = i->entity;

    const bool locally_destroyed = !is_valid(entity);
    const bool locally_owned_and_modified = locally_destroyed ? false :
      Modified == state(entity) &&
      m_parallel_rank   == i->owner ;

    if ( locally_destroyed || locally_owned_and_modified ) {

      // m_ghosting[0] is the SHARED communication

      for ( size_t j = ghosting_count ; --j ; ) {
        if ( entity_comm_map_erase( i->key, *m_ghosting[j] ) ) {
          ghosting_change_flags[ j ] = true ;
        }
      }
    }
  }

  std::vector< int > ghosting_change_flags_global( ghosting_count , 0 );

  BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "      internal_resolve_ghosted_modify_delete calling all_reduce_sum");
  all_reduce_sum( m_parallel_machine ,
                  & ghosting_change_flags[0] ,
                  & ghosting_change_flags_global[0] ,
                  ghosting_change_flags.size() );

  for ( unsigned ic = 0 ; ic < ghosting_change_flags_global.size() ; ++ic ) {
    if ( ghosting_change_flags_global[ic] ) {
      m_ghosting[ic]->m_sync_count = m_sync_count ;
    }
  }
  // BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "  exiting internal_resolve_ghosted_modify_delete");
}

void BulkData::resolve_ownership_of_modified_entities( const std::vector<Entity> &shared_modified )
{
    CommAll comm_all( m_parallel_machine );

    for ( int phase = 0; phase < 2; ++phase ) {
        for ( std::vector<Entity>::const_iterator i = shared_modified.begin() ; i != shared_modified.end() ; ++i ) {
            Entity entity = *i ;
            if ( parallel_owner_rank(entity) == m_parallel_rank &&
                   state(entity)  != Created ) {
                for ( PairIterEntityComm jc = entity_comm_map_shared(entity_key(entity)) ; ! jc.empty() ; ++jc ) {
                    comm_all.send_buffer( jc->proc ) .pack<EntityKey>( entity_key(entity) );
                }
            }
        }

        if (phase == 0) { //allocation phase
            BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "      calling allocate_buffers from internal_resolve_parallel_create");
            comm_all.allocate_buffers( m_parallel_size / 4 );
        }
        else { // communication phase
            BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "      calling communicate from internal_resolve_parallel_create");
            comm_all.communicate();
        }
    }

    for ( int receive_proc = 0 ; receive_proc < m_parallel_size ; ++receive_proc ) {
        CommBuffer & buf = comm_all.recv_buffer( receive_proc );
        EntityKey key ;
        while ( buf.remaining() ) {
            buf.unpack<EntityKey>( key );
            Entity entity = get_entity( key );

            // Set owner, will correct part membership later
            const bool changed = this->set_parallel_owner_rank_but_not_comm_lists( entity, receive_proc);
            if (changed) {
                internal_change_owner_in_comm_data(key, receive_proc);
            }
        }
    }
}

void BulkData::move_entities_to_proper_part_ownership( const std::vector<Entity> &shared_modified )
{
    std::ostringstream error_msg ;
    int error_flag = 0 ;

    PartVector shared_part , owned_part ;
    shared_part.push_back( & m_mesh_meta_data.globally_shared_part() );
    owned_part.push_back(  & m_mesh_meta_data.locally_owned_part() );

    std::vector<Entity>::const_reverse_iterator iend = shared_modified.rend();
    for ( std::vector<Entity>::const_reverse_iterator
          i = shared_modified.rbegin() ; i != iend ; ++i)
    {
      Entity entity = *i ;

      if ( parallel_owner_rank(entity) == m_parallel_rank &&
           state(entity) == Created )
      {
        // Created and not claimed by an existing owner

        const int new_owner = determine_new_owner( entity );

        const bool changed = this->set_parallel_owner_rank_but_not_comm_lists( entity, new_owner);
        if (changed) {
          internal_change_owner_in_comm_data(entity_key(entity), new_owner);
        }
      }

      if ( parallel_owner_rank(entity) != m_parallel_rank ) {
        // Do not own it and still have it.
        // Remove the locally owned, add the globally_shared
        set_synchronized_count( entity, m_sync_count);
        internal_change_entity_parts( entity , shared_part /*add*/, owned_part /*remove*/);
      }
      else if ( ! entity_comm_map_shared(entity_key(entity)).empty() ) {
        // Own it and has sharing information.
        // Add the globally_shared
        internal_change_entity_parts( entity , shared_part /*add*/, PartVector() /*remove*/ );
      }
      else {
        // Own it and does not have sharing information.
        // Remove the globally_shared
        internal_change_entity_parts( entity , PartVector() /*add*/, shared_part /*remove*/);
      }

      // Newly created shared entity had better be in the owned closure
      if ( ! in_owned_closure(*this, entity , m_parallel_rank ) ) {
        if ( 0 == error_flag ) {
          error_flag = 1 ;
          error_msg
            << "\nP" << m_parallel_rank << ": " << " FAILED\n"
            << "  The following entities were declared on multiple processors,\n"
            << "  cannot be parallel-shared, and were declared with"
            << "  parallel-ghosting information. {\n";
        }
        error_msg << "    " << print_entity_key(m_mesh_meta_data, entity_key(entity));
        error_msg << " also declared on" ;
        for ( PairIterEntityComm ec = entity_comm_map_shared(entity_key(entity)); ! ec.empty() ; ++ec ) {
          error_msg << " P" << ec->proc ;
        }
        error_msg << "\n" ;
      }
    }

    // Parallel-consistent error checking of above loop
    if ( error_flag ) { error_msg << "}\n" ; }
    all_reduce( m_parallel_machine , ReduceMax<1>( & error_flag ) );
    ThrowErrorMsgIf( error_flag, error_msg.str() );
}


void BulkData::update_comm_list(const std::vector<stk::mesh::Entity>& sharedModifiedEntities)
{
    // ------------------------------------------------------------
    // Update m_entity_comm based on shared_modified

    const size_t n_old = m_entity_comm_list.size();

    m_entity_comm_list.reserve(m_entity_comm_list.size() + sharedModifiedEntities.size());
    for (size_t i = 0, e = sharedModifiedEntities.size(); i < e; ++i) {
      Entity entity = sharedModifiedEntities[i];
      EntityCommListInfo new_comm = {entity_key(entity), entity, parallel_owner_rank(entity), NULL};
      m_entity_comm_list.push_back(new_comm);
    }

    std::sort( m_entity_comm_list.begin() + n_old , m_entity_comm_list.end() );

    std::inplace_merge( m_entity_comm_list.begin() ,
                        m_entity_comm_list.begin() + n_old ,
                        m_entity_comm_list.end() );

    EntityCommListInfoVector::iterator iter =
      std::unique( m_entity_comm_list.begin() , m_entity_comm_list.end() );

    m_entity_comm_list.erase( iter , m_entity_comm_list.end() );

    for(size_t i=0; i<m_entity_comm_list.size(); ++i) {
      EntityKey key = m_entity_comm_list[i].key;
      const EntityComm* entity_comm = m_entity_comm_map.entity_comm(key);
      m_entity_comm_list[i].entity_comm = entity_comm;
    }

    internal_sync_comm_list_owners();
}
//----------------------------------------------------------------------

// Postconditions:
//  * All shared entities have parallel-consistent owner
//  * Part membership of shared entities is up-to-date
//  * m_entity_comm is up-to-date
void BulkData::internal_resolve_parallel_create()
{
  Trace_("stk::mesh::BulkData::internal_resolve_parallel_create");
  BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "  entered internal_resolve_parallel_create");

  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");
  std::vector<Entity> shared_modified ;

  // Update the parallel index and
  // output shared and modified entities.
  internal_update_distributed_index(shared_modified );

  // ------------------------------------------------------------
  // Claim ownership on all shared_modified entities that I own
  // and which were not created in this modification cycle. All
  // sharing procs will need to be informed of this claim.

  resolve_ownership_of_modified_entities( shared_modified );

  // ------------------------------------------------------------
  // Update shared created entities.
  // - Revise ownership to selected processor
  // - Update sharing.
  // - Work backward so the 'in_owned_closure' function
  //   can evaluate related higher ranking entities.

  move_entities_to_proper_part_ownership( shared_modified );

  update_comm_list( shared_modified );
}



//----------------------------------------------------------------------

namespace {

#ifdef STK_VERBOSE_OUTPUT

bool no_buckets(const stk::mesh::BulkData& mesh, const stk::mesh::Part& part)
{
  for(stk::topology::rank_t r = stk::topology::NODE_RANK; r < mesh.mesh_meta_data().entity_rank_count(); ++r) {
    stk::mesh::Selector selector = part;
    const stk::mesh::BucketVector& buckets = mesh.get_buckets(r, selector);
    if (buckets.size() > 0) {
      return false;
    }
  }

  return true;
}

void print_bucket_data(const stk::mesh::BulkData& mesh)
{
  const stk::mesh::PartVector& all_parts = mesh.mesh_meta_data().get_parts();
  for(size_t i=0; i<all_parts.size(); ++i) {
    std::cout << "Part: " << all_parts[i]->name()<<std::endl;
    if (no_buckets(mesh, *all_parts[i])) {
      std::cout<<"\tEmpty"<<std::endl;
      continue;
    }
    for(stk::topology::rank_t r = stk::topology::NODE_RANK; r < mesh.mesh_meta_data().entity_rank_count(); ++r) {
      stk::mesh::Selector selector = *all_parts[i];
      const stk::mesh::BucketVector& buckets = mesh.get_buckets(r, selector);
      std::cout<<"\t"<< buckets.size() << " "<< r << " buckets";
      size_t min_entities = 1000000, max_entities = 0, total_entities = 0;
      double avg_entities = 0.0;
      for(size_t j=0; j<buckets.size(); ++j) {
        total_entities += buckets[j]->size();
        min_entities = std::min(min_entities, buckets[j]->size());
        max_entities = std::max(max_entities, buckets[j]->size());
      }
      avg_entities = buckets.size()>0 ? (1.0*total_entities)/buckets.size() : 0.0;
      if (total_entities == 0) {
        min_entities = 0;
      }
      std::cout << "; min=" << min_entities << ", avg=" << avg_entities << ", max=" << max_entities << ", tot=" << total_entities << std::endl;
    }
  }
}

#endif

}

//----------------------------------------------------------------------

bool BulkData::modification_end( modification_optimization opt)
{
  Trace_("stk::mesh::BulkData::modification_end");

  bool return_value = internal_modification_end( true, opt );

#ifdef GATHER_GET_BUCKETS_METRICS
  ++m_num_modifications;
#endif

#ifdef STK_VERBOSE_OUTPUT
  print_bucket_data(*this);
#endif

  m_modification_begin_description = "NO_MODIFICATION";

  write_modification_counts();

  return return_value;
}

bool BulkData::modification_end_for_entity_creation( EntityRank entity_rank, modification_optimization opt)
{
  Trace_("stk::mesh::BulkData::modification_end");

  bool return_value = internal_modification_end_for_entity_creation( entity_rank, true, opt );

#ifdef GATHER_GET_BUCKETS_METRICS
  ++m_num_modifications;
#endif

#ifdef STK_VERBOSE_OUTPUT
  print_bucket_data(*this);
#endif

  return return_value;
}

#if 0

namespace {

// Very, very handy for debugging parallel resolution...

void print_comm_list( const BulkData & mesh , bool doit )
{
  if ( doit ) {
    std::ostringstream msg ;

    msg << std::endl ;

    for ( EntityCommListInfoVector::const_iterator
          i =  mesh.comm_list().begin() ;
          i != mesh.comm_list().end() ; ++i ) {

      Entity entity = i->entity;
      msg << "P" << mesh.parallel_rank() << ": " ;

      print_entity_key( msg , MetaData::get(mesh) , i->key );

      msg << " owner(" << i->owner << ")" ;

      if ( !entity.is_valid() ) { msg << " del" ; }
      else if ( Modified == entity.state() ) { msg << " mod" ; }
      else { msg << "    " ; }

      for ( PairIterEntityComm ec = mesh.comm_list(i->key); ! ec.empty() ; ++ec ) {
        msg << " (" << ec->ghost_id << "," << ec->proc << ")" ;
      }
      msg << std::endl ;
    }

    std::cout << msg.str();
  }
}

}

#endif

void BulkData::update_comm_list_based_on_changes_in_comm_map()
// Resolution of shared and ghost modifications can empty
// the communication information for entities.
// If there is no communication information then the
// entity must be removed from the communication list.
{
  EntityCommListInfoVector::iterator i = m_entity_comm_list.begin();
  bool changed = false ;
  for ( ; i != m_entity_comm_list.end() ; ++i ) {
    const EntityComm* entity_comm = m_entity_comm_map.entity_comm(i->key);
    if ( entity_comm == NULL || entity_comm->comm_map.empty() ) {
      i->key = EntityKey();
      changed = true;
    }
    i->entity_comm = entity_comm;
  }
  if ( changed ) {
    i = std::remove_if( m_entity_comm_list.begin() ,
                        m_entity_comm_list.end() , IsInvalid() );
    m_entity_comm_list.erase( i , m_entity_comm_list.end() );
  }
}

bool BulkData::internal_modification_end_for_change_entity_owner( bool regenerate_aura, modification_optimization opt )
{
  Trace_("stk::mesh::BulkData::internal_modification_end");

  // The two states are MODIFIABLE and SYNCHRONiZED
  if ( m_sync_state == SYNCHRONIZED ) { return false ; }

  ThrowAssertMsg(impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

  ThrowAssertMsg(add_fmwk_data() || impl::check_no_shared_elements_or_higher(*this)==0, "BulkData::modification_end ERROR, Sharing of entities with rank ELEMENT_RANK or higher is not allowed.");

  if (parallel_size() > 1) {
    // Resolve modification or deletion of shared entities
    // which can cause deletion of ghost entities.

    //passes np1 segfaults (np2-np4) if not called
    internal_resolve_shared_modify_delete();

    // Resolve modification or deletion of ghost entities
    // by destroying ghost entities that have been touched.

    //doesn't fail if not called (np1-np4)
    //internal_resolve_ghosted_modify_delete();

    //passes np1 segfaults np2-np4 if not called
    update_comm_list_based_on_changes_in_comm_map();

    // Resolve creation of entities: discover sharing and set unique ownership.
    //passes np1 fails on 10 tests np2, 3 tests np3, segfaults np4
    internal_resolve_parallel_create();

    // Resolve part membership for shared entities.
    // This occurs after resolving creation so created and shared
    // entities are resolved along with previously existing shared entities.

    //doesn't fail if not called (np1-np4)
    //internal_resolve_shared_membership();

    // Regenerate the ghosting aura around all shared mesh entities.
    if ( regenerate_aura )
    {
      //passes on np1 hangs on np2-np4 if not called
      internal_regenerate_aura();
    }

    // ------------------------------
    // Verify parallel consistency of mesh entities.
    // Unique ownership, communication lists, sharing part membership,
    // application part membership consistency.
#ifndef NDEBUG
    std::ostringstream msg ;
    bool is_consistent = true;
    is_consistent = impl::comm_mesh_verify_parallel_consistency( *this , msg );
    ThrowErrorMsgIf( !is_consistent, msg.str() );
#endif
  }
  else {
      std::vector<Entity> shared_modified ;
      //doesn't fail if not called (np1-np4)
      //internal_update_distributed_index( shared_modified );
  }

  // ------------------------------
  // Now sort the bucket entities.
  // This does not change the entities, relations, or field data.
  // However, it insures that the ordering of entities and buckets
  // is independent of the order in which a set of changes were
  // performed.
  //
  //optimize_buckets combines multiple buckets in a bucket-family into
  //a single larger bucket, and also does a sort.
  //If optimize_buckets has not been requested, still do the sort.

  if ( opt == MOD_END_COMPRESS_AND_SORT ) {
    //doesn't fail if not called (np1-np4)
    //m_bucket_repository.optimize_buckets();
  }
  else {
    //doesn't fail if not called (np1-np4)
    //m_bucket_repository.internal_sort_bucket_entities();
  }

  // ------------------------------

  //doesn't fail if not called (np1-np4)
  //m_bucket_repository.internal_modification_end();

  //doesn't fail if not called (np1-np4)
  //internal_update_fast_comm_maps();

  m_sync_state = SYNCHRONIZED ;
  m_add_node_sharing_called = false;

  //doesn't fail if not called (np1-np4)
  //update_deleted_entities_container();

  return true ;
}

bool BulkData::internal_modification_end( bool regenerate_aura, modification_optimization opt )
{
  Trace_("stk::mesh::BulkData::internal_modification_end");

  // The two states are MODIFIABLE and SYNCHRONiZED
  if ( m_sync_state == SYNCHRONIZED ) { return false ; }

  ThrowAssertMsg(impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

  ThrowAssertMsg(add_fmwk_data() || impl::check_no_shared_elements_or_higher(*this)==0, "BulkData::modification_end ERROR, Sharing of entities with rank ELEMENT_RANK or higher is not allowed.");

  if (parallel_size() > 1) {
    // Resolve modification or deletion of shared entities
    // which can cause deletion of ghost entities.
    internal_resolve_shared_modify_delete();

    // Resolve modification or deletion of ghost entities
    // by destroying ghost entities that have been touched.
    internal_resolve_ghosted_modify_delete();

    update_comm_list_based_on_changes_in_comm_map();

    // Resolve creation of entities: discover sharing and set unique ownership.
    internal_resolve_parallel_create();

    // Resolve part membership for shared entities.
    // This occurs after resolving creation so created and shared
    // entities are resolved along with previously existing shared entities.
    internal_resolve_shared_membership();

    // Regenerate the ghosting aura around all shared mesh entities.
    if ( regenerate_aura ) { internal_regenerate_aura(); }

    // ------------------------------
    // Verify parallel consistency of mesh entities.
    // Unique ownership, communication lists, sharing part membership,
    // application part membership consistency.
#ifndef NDEBUG
    std::ostringstream msg ;
    bool is_consistent = true;
    is_consistent = impl::comm_mesh_verify_parallel_consistency( *this , msg );
    ThrowErrorMsgIf( !is_consistent, msg.str() );
#endif
  }
  else {
      std::vector<Entity> shared_modified ;
      internal_update_distributed_index( shared_modified );
  }

  // ------------------------------
  // Now sort the bucket entities.
  // This does not change the entities, relations, or field data.
  // However, it insures that the ordering of entities and buckets
  // is independent of the order in which a set of changes were
  // performed.
  //
  //optimize_buckets combines multiple buckets in a bucket-family into
  //a single larger bucket, and also does a sort.
  //If optimize_buckets has not been requested, still do the sort.

  if ( opt == MOD_END_COMPRESS_AND_SORT ) {
    m_bucket_repository.optimize_buckets();
  }
  else {
    m_bucket_repository.internal_sort_bucket_entities();
  }

  // ------------------------------

  m_bucket_repository.internal_modification_end();

  internal_update_fast_comm_maps();

  m_sync_state = SYNCHRONIZED ;
  m_add_node_sharing_called = false;

  update_deleted_entities_container();

  return true ;
}

//////////////////////////////////// Free funcions to help with modification end (exp) for edges

void fillElementsConnectedToNodes(stk::mesh::BulkData &stkMeshBulkData, std::vector<stk::mesh::Entity> &nodes,
                                  std::vector<stk::mesh::Entity> & elementsConnectedToNodes)
{
    stk::mesh::Entity const * elemStartNode = stkMeshBulkData.begin_elements(nodes[0]);
    stk::mesh::Entity const * elemEndNode = stkMeshBulkData.end_elements(nodes[0]);
    elementsConnectedToNodes.assign(elemStartNode, elemEndNode);
    std::sort(elementsConnectedToNodes.begin(), elementsConnectedToNodes.end());
    std::vector<stk::mesh::Entity>::iterator intersectionEnd = elementsConnectedToNodes.end();

    std::vector<stk::mesh::Entity> elementsConnectedToNodesTemp;
    for (size_t i = 1; i < nodes.size(); ++i)
    {
      elementsConnectedToNodesTemp.clear();
      elementsConnectedToNodesTemp.assign(elementsConnectedToNodes.begin(), intersectionEnd);
      elemStartNode = stkMeshBulkData.begin_elements(nodes[i]);
      elemEndNode = stkMeshBulkData.end_elements(nodes[i]);
      std::vector<stk::mesh::Entity> elems(elemStartNode, elemEndNode);
      std::sort(elems.begin(), elems.end());

      intersectionEnd = std::set_intersection( elems.begin(), elems.end(),
                                               elementsConnectedToNodesTemp.begin(), elementsConnectedToNodesTemp.end(),
                                               elementsConnectedToNodes.begin() );
      if (intersectionEnd == elementsConnectedToNodes.begin()) break;  // Empty set
    }

    elementsConnectedToNodes.resize(intersectionEnd-elementsConnectedToNodes.begin());
}

void fillFacesConnectedToNodes(stk::mesh::BulkData &stkMeshBulkData, std::vector<stk::mesh::Entity> &nodes,
              std::vector<stk::mesh::Entity> & facesConnectedToNodes)
{
    facesConnectedToNodes.clear();
    stk::mesh::EntityRank needConnectivityOfType = stk::topology::FACE_RANK;
    if (stkMeshBulkData.connectivity_map().valid(stk::topology::NODE_RANK, needConnectivityOfType))
    {
        //TODO: fill faces for all edge nodes
        if ( stkMeshBulkData.num_faces(nodes[0]) > 0 && stkMeshBulkData.num_faces(nodes[1]) > 0 )
        {
            facesConnectedToNodes.resize(20);
            stk::mesh::Entity const * facesStartNode1 = stkMeshBulkData.begin_faces(nodes[0]);
            stk::mesh::Entity const * facesEndNode1 = stkMeshBulkData.end_faces(nodes[0]);
            stk::mesh::Entity const * facesStartNode2 = stkMeshBulkData.begin_faces(nodes[1]);
            stk::mesh::Entity const * facesEndNode2 = stkMeshBulkData.end_faces(nodes[1]);

            std::vector<stk::mesh::Entity> faces1(facesStartNode1, facesEndNode1);
            std::sort(faces1.begin(), faces1.end());
            std::vector<stk::mesh::Entity> faces2(facesStartNode2, facesEndNode2);
            std::sort(faces2.begin(), faces2.end());

            std::vector<stk::mesh::Entity>::iterator iter = std::set_intersection( faces1.begin(), faces1.end(),
                  faces2.begin(), faces2.end(), facesConnectedToNodes.begin());

            facesConnectedToNodes.resize(iter-facesConnectedToNodes.begin());
        }
    }
    else
    {
        stk::mesh::EntityVector entitiesNode1;
        stk::mesh::EntityVector entitiesNode2;
        int numFacesConnectedToNode1 = stk::mesh::get_connectivity(stkMeshBulkData, nodes[0], needConnectivityOfType, entitiesNode1);
        int numFacesConnectedToNode2 = stk::mesh::get_connectivity(stkMeshBulkData, nodes[1], needConnectivityOfType, entitiesNode2);
        if ( numFacesConnectedToNode1 > 0 && numFacesConnectedToNode2 > 0 )
        {
            facesConnectedToNodes.resize(entitiesNode1.size()+entitiesNode2.size());

            std::vector<stk::mesh::Entity>::iterator iter = std::set_intersection( entitiesNode1.begin(), entitiesNode1.end(),
                    entitiesNode2.begin(), entitiesNode2.end(), facesConnectedToNodes.begin());

            facesConnectedToNodes.resize(iter-facesConnectedToNodes.begin());
        }
    }
}

bool doesEdgeNeedGhostingCommunication(stk::mesh::BulkData &stkMeshBulkData, std::vector<stk::mesh::Entity>& connectedEntities)
{
    bool communicate_edge_for_ghosting = false;
    for (size_t j=0;j<connectedEntities.size();j++)
    {
        bool isEntityGhostedToAnotherProc = stkMeshBulkData.is_ghosted_onto_another_proc(stkMeshBulkData.entity_key(connectedEntities[j]));
        if ( isEntityGhostedToAnotherProc )
        {
            communicate_edge_for_ghosting = true;
        }
        else
        {
            connectedEntities[j] = Entity();
        }
    }
    return communicate_edge_for_ghosting;
}


void connectUpwardEntityToEntity(stk::mesh::BulkData& mesh, stk::mesh::Entity upward_entity,
        stk::mesh::Entity entity, std::vector<stk::mesh::Entity> &nodes)
{
    // get node upwardentity ids
    uint num_nodes = mesh.num_nodes(entity);
    EntityRank entity_rank = mesh.entity_rank(entity);
    std::vector<stk::mesh::EntityId> nodeIdsForEntity(num_nodes);
    for (uint i = 0 ; i < num_nodes ; ++i)
    {
      nodeIdsForEntity[i]= mesh.identifier(nodes[i]);
    }

    // scratch space
    stk::mesh::OrdinalVector ordinal_scratch;
    ordinal_scratch.reserve(64);
    stk::mesh::PartVector part_scratch;
    part_scratch.reserve(64);
    stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(0);

    // now what
    stk::topology upward_entity_topology = mesh.bucket(upward_entity).topology();
    std::vector<stk::mesh::EntityId> nodeIds(num_nodes);
    std::vector<stk::mesh::EntityId> upward_entity_node_ids(upward_entity_topology.num_nodes());
    unsigned entity_ordinal = 100000;
    stk::mesh::Entity const * upward_entity_nodes = mesh.begin_nodes(upward_entity);
    for (size_t k=0;k<upward_entity_topology.num_nodes();k++)
    {
        upward_entity_node_ids[k] = mesh.identifier(upward_entity_nodes[k]);
    }

    stk::topology entity_top;
    for (size_t k=0;k<upward_entity_topology.num_sub_topology(entity_rank);k++)
    {
        if(entity_rank == stk::topology::EDGE_RANK)
        {
          upward_entity_topology.edge_nodes(upward_entity_node_ids, k, nodeIds.begin());
          entity_top = upward_entity_topology.edge_topology();
        }
        else
        {
          upward_entity_topology.face_nodes(upward_entity_node_ids, k, nodeIds.begin());
          entity_top = upward_entity_topology.face_topology();
        }
        if ( entity_top.equivalent(nodeIds, nodeIdsForEntity).first )
        {
            entity_ordinal = k;
            break;
        }
    }
    ThrowRequireMsg(entity_ordinal !=100000, "Program error. Contact sierra-help for support.");
    mesh.declare_relation(upward_entity, entity, entity_ordinal, perm, ordinal_scratch, part_scratch);
}

void connectGhostedEntitiesToEdge(stk::mesh::BulkData &stkMeshBulkData, std::vector<stk::mesh::Entity> &entitiesConnectedToNodes, stk::mesh::Entity edge, std::vector<stk::mesh::Entity>& nodes)
{
    for (size_t j=0; j<entitiesConnectedToNodes.size();j++)
    {
        bool isEntityGhostedOntoThisProc = stkMeshBulkData.in_receive_ghost(stkMeshBulkData.aura_ghosting(), stkMeshBulkData.entity_key(entitiesConnectedToNodes[j]));

        if ( isEntityGhostedOntoThisProc )
        {
            impl::connectEntityToEdge(stkMeshBulkData, entitiesConnectedToNodes[j], edge, nodes);
        }
    }
}

void connectGhostedEntitiesToEntity(stk::mesh::BulkData &stkMeshBulkData, std::vector<stk::mesh::Entity> &entitiesConnectedToNodes, stk::mesh::Entity entity, std::vector<stk::mesh::Entity>& nodes)
{
    for (size_t j=0; j<entitiesConnectedToNodes.size();j++)
    {
        bool isEntityGhostedOntoThisProc = stkMeshBulkData.in_receive_ghost(stkMeshBulkData.aura_ghosting(), stkMeshBulkData.entity_key(entitiesConnectedToNodes[j]));

        if ( isEntityGhostedOntoThisProc )
        {
            connectUpwardEntityToEntity(stkMeshBulkData, entitiesConnectedToNodes[j], entity, nodes);
        }
    }
}

void determineEntitiesThatNeedGhosting(stk::mesh::BulkData &stkMeshBulkData, stk::mesh::Entity edge, std::vector<stk::mesh::Entity>& entitiesConnectedToNodes,
        std::vector<stk::mesh::Entity> &nodes, std::set<EntityProc, EntityLess> &addGhostedEntities)
{
    // Grab all the entities attached to the 2 nodes
    // If the entity is ghosted and the edge is owned, then the edge needs to be ghosted.
    bool doesEdgeNeedToBeGhosted = doesEdgeNeedGhostingCommunication(stkMeshBulkData, entitiesConnectedToNodes);
    if ( doesEdgeNeedToBeGhosted )
    {
        for (size_t j=0;j<entitiesConnectedToNodes.size();j++)
        {
            if ( entitiesConnectedToNodes[j] != Entity() )
            {
                PairIterEntityComm ghosted = stkMeshBulkData.entity_comm_map( stkMeshBulkData.entity_key(entitiesConnectedToNodes[j]) , stkMeshBulkData.aura_ghosting());
                for (PairIterEntityComm ec = ghosted; !ec.empty(); ++ec)
                {
                    if ( ec->proc != stkMeshBulkData.parallel_rank() )
                    {
                        bool isEdgeSharedWithOtherProc = stkMeshBulkData.in_shared(stkMeshBulkData.entity_key(edge), ec->proc);
                        bool isEntityGhostedToAnotherProc = stkMeshBulkData.is_ghosted_onto_proc(stkMeshBulkData.entity_key(entitiesConnectedToNodes[j]), ec->proc);

                        if ( !isEdgeSharedWithOtherProc && isEntityGhostedToAnotherProc )
                        {
                            addGhostedEntities.insert(EntityProc(edge, ec->proc));
                        }
                    }
                }
            }
        }
    }
}

void connect_edges_to_ghosted_entities_and_fill_edges_to_ghost_onto_other_processors(stk::mesh::BulkData &stkMeshBulkata, std::set<EntityProc, EntityLess> &entitiesToGhostOntoOtherProcessors)
{
    const stk::mesh::BucketVector& edge_buckets = stkMeshBulkata.buckets(stk::topology::EDGE_RANK);

    for(size_t bucketIndex = 0; bucketIndex < edge_buckets.size(); bucketIndex++)
    {
        const stk::mesh::Bucket& bucket = *edge_buckets[bucketIndex];
        for(size_t entityIndex = 0; entityIndex < bucket.size(); entityIndex++)
        {
            Entity edge = bucket[entityIndex];
            if ( stkMeshBulkata.state(edge) == Unchanged ) continue;

            std::vector<stk::mesh::Entity> edge_nodes(2);
            stk::mesh::Entity const *  nodes = stkMeshBulkata.begin_nodes(edge);
            edge_nodes[0] = nodes[0];
            edge_nodes[1] = nodes[1];

            std::vector<stk::mesh::Entity> facesConnectedToNodes;
            fillFacesConnectedToNodes(stkMeshBulkata, edge_nodes, facesConnectedToNodes);
            connectGhostedEntitiesToEdge(stkMeshBulkata, facesConnectedToNodes, edge, edge_nodes);

            std::vector<stk::mesh::Entity> elementsConnectedToNodes;
            fillElementsConnectedToNodes(stkMeshBulkata, edge_nodes, elementsConnectedToNodes);
            connectGhostedEntitiesToEdge(stkMeshBulkata, elementsConnectedToNodes, edge, edge_nodes);

            if ( bucket.owned() || bucket.shared() )
            {
                determineEntitiesThatNeedGhosting(stkMeshBulkata, edge, facesConnectedToNodes, edge_nodes, entitiesToGhostOntoOtherProcessors);
                determineEntitiesThatNeedGhosting(stkMeshBulkata, edge, elementsConnectedToNodes, edge_nodes, entitiesToGhostOntoOtherProcessors);
            }
        }
    }

    std::set< EntityKey > entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs;
    comm_sync_send_recv(stkMeshBulkata, entitiesToGhostOntoOtherProcessors, entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs);
}

void find_upward_connected_entities_to_ghost_onto_other_processors(stk::mesh::BulkData &mesh, std::set<EntityProc, EntityLess> &entitiesToGhostOntoOtherProcessors, EntityRank entity_rank)
{
    const stk::mesh::BucketVector& entity_buckets = mesh.buckets(entity_rank);
    bool isedge = (entity_rank == stk::topology::EDGE_RANK);

    for(size_t bucketIndex = 0; bucketIndex < entity_buckets.size(); bucketIndex++)
    {
        const stk::mesh::Bucket& bucket = *entity_buckets[bucketIndex];
        size_t num_nodes_on_entity = bucket.topology().num_nodes();
        std::vector<stk::mesh::Entity> entity_nodes(num_nodes_on_entity);

        for(size_t entityIndex = 0; entityIndex < bucket.size(); entityIndex++)
        {
            Entity entity = bucket[entityIndex];
            if ( mesh.state(entity) == Unchanged ) continue;

            entity_nodes.assign(mesh.begin_nodes(entity), mesh.end_nodes(entity));
            std::vector<stk::mesh::Entity> facesConnectedToNodes;

            if(isedge)
            {
              fillFacesConnectedToNodes(mesh, entity_nodes, facesConnectedToNodes);
              //connectGhostedEntitiesToEdge(mesh, facesConnectedToNodes, entity, entity_nodes);
              connectGhostedEntitiesToEntity(mesh, facesConnectedToNodes, entity, entity_nodes);
            }

            std::vector<stk::mesh::Entity> elementsConnectedToNodes;
            fillElementsConnectedToNodes(mesh, entity_nodes, elementsConnectedToNodes);
            //connectGhostedEntitiesToEdge(mesh, elementsConnectedToNodes, entity, entity_nodes);
            connectGhostedEntitiesToEntity(mesh, elementsConnectedToNodes, entity, entity_nodes);

            if ( bucket.owned() || bucket.shared() )
            {
                if (isedge)
                {
                  determineEntitiesThatNeedGhosting(mesh, entity, facesConnectedToNodes, entity_nodes, entitiesToGhostOntoOtherProcessors);
                }
                determineEntitiesThatNeedGhosting(mesh, entity, elementsConnectedToNodes, entity_nodes, entitiesToGhostOntoOtherProcessors);
            }
        }
    }

    std::set< EntityKey > entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs;
    comm_sync_send_recv(mesh, entitiesToGhostOntoOtherProcessors, entitiesGhostedOnThisProcThatNeedInfoFromOtherProcs);
}

void connect_ghosted_edges_received_to_ghosted_faces_and_elements(stk::mesh::BulkData &stkMeshBulkData)
{
    if (!stkMeshBulkData.connectivity_map().valid(stk::topology::EDGE_RANK, stk::topology::FACE_RANK))
    {
        const stk::mesh::BucketVector& edge_buckets = stkMeshBulkData.buckets(stk::topology::EDGE_RANK);
        for(size_t bucketIndex = 0; bucketIndex < edge_buckets.size(); bucketIndex++)
        {
            const stk::mesh::Bucket& bucket = *edge_buckets[bucketIndex];
            if ( bucket.in_aura() )
            {
                for(size_t entityIndex = 0; entityIndex < bucket.size(); entityIndex++)
                {
                    Entity edge = bucket[entityIndex];
                    if ( stkMeshBulkData.state(edge) == Unchanged ) continue;

                    std::vector<stk::mesh::Entity> edge_nodes(2);
                    stk::mesh::Entity const *  nodes = stkMeshBulkData.begin_nodes(edge);
                    edge_nodes[0] = nodes[0];
                    edge_nodes[1] = nodes[1];

                    std::vector<stk::mesh::Entity> facesConnectedToNodes;
                    fillFacesConnectedToNodes(stkMeshBulkData, edge_nodes, facesConnectedToNodes);
                    connectGhostedEntitiesToEdge(stkMeshBulkData, facesConnectedToNodes, edge, edge_nodes);
                }
            }
        }
    }

    if (!stkMeshBulkData.connectivity_map().valid(stk::topology::EDGE_RANK, stk::topology::ELEMENT_RANK) )
    {
        const stk::mesh::BucketVector& edge_buckets = stkMeshBulkData.buckets(stk::topology::EDGE_RANK);
        for(size_t bucketIndex = 0; bucketIndex < edge_buckets.size(); bucketIndex++)
        {
            const stk::mesh::Bucket& bucket = *edge_buckets[bucketIndex];
            if ( bucket.in_aura() )
            {
                for(size_t entityIndex = 0; entityIndex < bucket.size(); entityIndex++)
                {
                    Entity edge = bucket[entityIndex];
                    if ( stkMeshBulkData.state(edge) == Unchanged ) continue;

                    std::vector<stk::mesh::Entity> edge_nodes(2);
                    stk::mesh::Entity const *  nodes = stkMeshBulkData.begin_nodes(edge);
                    edge_nodes[0] = nodes[0];
                    edge_nodes[1] = nodes[1];

                    std::vector<stk::mesh::Entity> elementsConnectedToNodes;
                    fillElementsConnectedToNodes(stkMeshBulkData, edge_nodes, elementsConnectedToNodes);
                    connectGhostedEntitiesToEdge(stkMeshBulkData, elementsConnectedToNodes, edge, edge_nodes);
                }
            }
        }
    }
}

void connect_ghosted_entities_received_to_ghosted_upwardly_connected_entities(stk::mesh::BulkData &mesh, EntityRank entity_rank)
{
    if (entity_rank == stk::topology::EDGE_RANK && !mesh.connectivity_map().valid(stk::topology::EDGE_RANK, stk::topology::FACE_RANK) )
    {
        const stk::mesh::BucketVector& entity_buckets = mesh.buckets(stk::topology::EDGE_RANK);
        for(size_t bucketIndex = 0; bucketIndex < entity_buckets.size(); bucketIndex++)
        {
            const stk::mesh::Bucket& bucket = *entity_buckets[bucketIndex];
            if ( bucket.in_aura() )
            {
                for(size_t entityIndex = 0; entityIndex < bucket.size(); entityIndex++)
                {
                    Entity entity = bucket[entityIndex];
                    if ( mesh.state(entity) == Unchanged ) continue;

                    std::vector<stk::mesh::Entity> entity_nodes(2);
                    stk::mesh::Entity const *  nodes = mesh.begin_nodes(entity);
                    entity_nodes[0] = nodes[0];
                    entity_nodes[1] = nodes[1];

                    std::vector<stk::mesh::Entity> facesConnectedToNodes;
                    fillFacesConnectedToNodes(mesh, entity_nodes, facesConnectedToNodes);
                    connectGhostedEntitiesToEdge(mesh, facesConnectedToNodes, entity, entity_nodes);
                }
            }
        }
    }

    if (!mesh.connectivity_map().valid(entity_rank, stk::topology::ELEMENT_RANK) )
    {
        const stk::mesh::BucketVector& entity_buckets = mesh.buckets(entity_rank);
        for(size_t bucketIndex = 0; bucketIndex < entity_buckets.size(); bucketIndex++)
        {
            const stk::mesh::Bucket& bucket = *entity_buckets[bucketIndex];
            if ( bucket.in_aura() )
            {
              size_t num_nodes_on_entity = bucket.topology().num_nodes();
              std::vector<stk::mesh::Entity> entity_nodes(num_nodes_on_entity);

                for(size_t entityIndex = 0; entityIndex < bucket.size(); entityIndex++)
                {
                    Entity entity = bucket[entityIndex];
                    if ( mesh.state(entity) == Unchanged ) continue;

                    entity_nodes.assign(mesh.begin_nodes(entity), mesh.end_nodes(entity));

                    std::vector<stk::mesh::Entity> elementsConnectedToNodes;
                    fillElementsConnectedToNodes(mesh, entity_nodes, elementsConnectedToNodes);
                    connectGhostedEntitiesToEdge(mesh, elementsConnectedToNodes, entity, entity_nodes);
                }
            }
        }
    }
}

bool BulkData::internal_modification_end_for_entity_creation( EntityRank entity_rank, bool regenerate_aura, modification_optimization opt )
{
  Trace_("stk::mesh::BulkData::internal_modification_end_for_entity_creation");

  // The two states are MODIFIABLE and SYNCHRONiZED
  if ( m_sync_state == SYNCHRONIZED ) { return false ; }

  ThrowAssertMsg(impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

  if (parallel_size() > 1)
  {
    std::vector<Entity> shared_modified ;

    // Update the parallel index and
    // output shared and modified entities.
    internal_update_distributed_index( entity_rank, shared_modified );

    // ------------------------------------------------------------
    // Claim ownership on all shared_modified entities that I own
    // and which were not created in this modification cycle. All
    // sharing procs will need to be informed of this claim.

    resolve_ownership_of_modified_entities( shared_modified );

    // ------------------------------------------------------------
    // Update shared created entities.
    // - Revise ownership to selected processor
    // - Update sharing.
    // - Work backward so the 'in_owned_closure' function
    //   can evaluate related higher ranking entities.

    move_entities_to_proper_part_ownership( shared_modified );

    update_comm_list( shared_modified );

    internal_resolve_shared_membership();

    std::set<EntityProc, EntityLess> entitiesToGhostOntoOtherProcessors(EntityLess(*this));

    //connect_edges_to_ghosted_entities_and_fill_edges_to_ghost_onto_other_processors(*this, entitiesToGhostOntoOtherProcessors);
    find_upward_connected_entities_to_ghost_onto_other_processors(*this, entitiesToGhostOntoOtherProcessors, entity_rank);

    ghost_entities_and_fields(aura_ghosting(), entitiesToGhostOntoOtherProcessors);

    //connect_ghosted_edges_received_to_ghosted_faces_and_elements(*this);
    connect_ghosted_entities_received_to_ghosted_upwardly_connected_entities(*this, entity_rank);

#ifndef NDEBUG
    std::ostringstream msg ;
    bool is_consistent = true;
    is_consistent = impl::comm_mesh_verify_parallel_consistency( *this , msg );
    ThrowErrorMsgIf( !is_consistent, msg.str() );
#endif
  }
  else {
      std::vector<Entity> shared_modified ;
      internal_update_distributed_index( entity_rank, shared_modified );
  }

  // ------------------------------
  // Now sort the bucket entities.
  // This does not change the entities, relations, or field data.
  // However, it insures that the ordering of entities and buckets
  // is independent of the order in which a set of changes were
  // performed.
  //
  //optimize_buckets combines multiple buckets in a bucket-family into
  //a single larger bucket, and also does a sort.
  //If optimize_buckets has not been requested, still do the sort.

  if ( opt == MOD_END_COMPRESS_AND_SORT ) {
    m_bucket_repository.optimize_buckets();
  }
  else {
    m_bucket_repository.internal_sort_bucket_entities();
  }

  // ------------------------------

  m_bucket_repository.internal_modification_end();

  internal_update_fast_comm_maps();

  m_sync_state = SYNCHRONIZED ;
  m_add_node_sharing_called = false;

  update_deleted_entities_container();

  return true ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

template <typename T>
T const* get_begin_itr(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank);

template <typename T>
T const* get_end_itr(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank);


template <>
Entity const* get_begin_itr<Entity>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin(bucket_ordinal, rank); }

template <>
Entity const* get_end_itr<Entity>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end(bucket_ordinal, rank); }


template <>
ConnectivityOrdinal const* get_begin_itr<ConnectivityOrdinal>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin_ordinals(bucket_ordinal, rank); }

template <>
ConnectivityOrdinal const* get_end_itr<ConnectivityOrdinal>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end_ordinals(bucket_ordinal, rank); }


template <>
Permutation const* get_begin_itr<Permutation>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin_permutations(bucket_ordinal, rank); }

template <>
Permutation const* get_end_itr<Permutation>(const Bucket& bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end_permutations(bucket_ordinal, rank); }

//
// Because the connectivity API in Bucket is not templated...
//

template <typename T>
T const *get_begin_relation_data(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank);

template <typename T>
T const *get_end_relation_data(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank);

template <>
Entity const *get_begin_relation_data<Entity>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin(bucket_ordinal, rank); }

template <>
Entity const *get_end_relation_data<Entity>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end(bucket_ordinal, rank); }

template <>
ConnectivityOrdinal const *get_begin_relation_data<ConnectivityOrdinal>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin_ordinals(bucket_ordinal, rank); }

template <>
ConnectivityOrdinal const *get_end_relation_data<ConnectivityOrdinal>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end_ordinals(bucket_ordinal, rank); }

template <>
Permutation const *get_begin_relation_data<Permutation>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.begin_permutations(bucket_ordinal, rank); }

template <>
Permutation const *get_end_relation_data<Permutation>(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{ return bucket.end_permutations(bucket_ordinal, rank); }

template <typename T>
void verify_relation_data(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank)
{
  T const * old_it = get_begin_relation_data<T>(bucket, bucket_ordinal, rank);
  T const * old_it_end = get_end_relation_data<T>(bucket, bucket_ordinal, rank);

  T const * new_data_begin = get_begin_itr<T>(bucket, bucket_ordinal, rank);
  T const * new_data_end   = get_end_itr<T>(bucket, bucket_ordinal, rank);

  ThrowRequire( std::distance(new_data_begin, new_data_end ) ==
                std::distance(old_it, old_it_end ) );

  ThrowRequire( std::distance(new_data_begin, new_data_end) ==
                bucket.num_connectivity(bucket_ordinal, rank) );

  T const * new_it = new_data_begin;
  for (; old_it != old_it_end ; ++old_it , ++new_it) {
    T old_data = *old_it;
    T new_data = *new_it;
    ThrowRequire(old_data == new_data);
  }
}

}

void BulkData::verify_relations(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank) const
{
  verify_relation_data<Entity>(bucket, bucket_ordinal, rank);
  verify_relation_data<ConnectivityOrdinal>(bucket, bucket_ordinal, rank);
  if (bucket.has_permutation(rank)) {
    verify_relation_data<Permutation>(bucket, bucket_ordinal, rank);
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

enum { PART_ORD_UNIVERSAL = 0 };
enum { PART_ORD_OWNED     = 1 };
enum { PART_ORD_SHARED    = 2 };

namespace {

void pack_induced_memberships( BulkData& bulk_data,
                               CommAll & comm ,
                               const EntityCommListInfoVector & entity_comm )
{
  OrdinalVector empty , induced ;
  for ( EntityCommListInfoVector::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    if ( bulk_data.in_shared( i->key , i->owner ) ) {
      // Is shared with owner, send to owner.

      empty.clear();
      induced.clear();

      induced_part_membership(bulk_data, i->entity , empty , induced );

      CommBuffer & buf = comm.send_buffer( i->owner );

      unsigned tmp = induced.size();

      buf.pack<unsigned>( tmp );

      for ( OrdinalVector::iterator
            j = induced.begin() ; j != induced.end() ; ++j ) {
        buf.pack<unsigned>( *j );
      }
    }
  }
}

void generate_send_list( BulkData& bulk_data,
                         const size_t sync_count ,
                         const int p_rank ,
                         const EntityCommListInfoVector & entity_comm ,
                               std::vector<EntityProc> & send_list )
{
  for ( EntityCommListInfoVector::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    if ( i->owner == p_rank &&
         bulk_data.synchronized_count(i->entity) == sync_count ) {

      for ( PairIterEntityComm ec = bulk_data.entity_comm_map(i->key); ! ec.empty(); ++ec ) {
        EntityProc tmp( i->entity , ec->proc );
        send_list.push_back( tmp );
      }
    }
  }

  {
    std::sort( send_list.begin() , send_list.end() , EntityLess(bulk_data) );
    std::vector<EntityProc>::iterator i =
      std::unique( send_list.begin() , send_list.end() );
    send_list.erase( i , send_list.end() );
  }
}

void pack_part_memberships( BulkData& meshbulk, CommAll & comm ,
                            const std::vector<EntityProc> & send_list )
{
  for ( std::vector<EntityProc>::const_iterator
        i = send_list.begin() ; i != send_list.end() ; ++i ) {

    Entity entity = i->first;

    std::pair<const unsigned *, const unsigned *>
      part_ord = meshbulk.bucket(entity).superset_part_ordinals();

    // I am the owner; therefore, the first three members are
    // universal, uses, and owns.  Don't send them.

    // I am the owner.  The first two memberships are
    // universal_part and locally_owned_part.  The third
    // membership may be globally_shared_part ;

    const unsigned count_all  = part_ord.second - part_ord.first ;
    const unsigned count_skip =
      ( 2 < count_all && part_ord.first[2] == PART_ORD_SHARED ) ? 3 : 2 ;

    const unsigned count_send = count_all - count_skip ;

    const unsigned * const start_send = part_ord.first + count_skip ;

    comm.send_buffer( i->second ).pack<EntityKey>( meshbulk.entity_key(entity) )
                                 .pack<unsigned>( count_send )
                                 .pack<unsigned>( start_send , count_send );
  }
}

}

//  Mesh entity membership changes must be synchronized among
//  processes that share mesh entities and propagated to
//  processes that ghost copies of the mesh entities.
//
//  Precondition: correct shared and ghosting lists.
//
//  Part memberships may have been added or removed
//  either explicitly or indirectly via entity relationships
//  being added or removed.q

void BulkData::internal_resolve_shared_membership()
{
    Trace_("stk::mesh::BulkData::internal_resolve_shared_membership");
    BABBLE_STK_PARALLEL_COMM(m_parallel_machine, "  entered internal_resolve_shared_membership");

    ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

    const MetaData & meta = m_mesh_meta_data;
    ParallelMachine p_comm = m_parallel_machine;
    const int p_rank = m_parallel_rank;
    const int p_size = m_parallel_size;
    const PartVector & all_parts = meta.get_parts();

    const Part & part_universal = meta.universal_part();
    const Part & part_owned = meta.locally_owned_part();
    const Part & part_shared = meta.globally_shared_part();

    // Quick verification of part ordinal assumptions

    ThrowRequireMsg(PART_ORD_UNIVERSAL == part_universal.mesh_meta_data_ordinal(),
            "Universal part ordinal is wrong, expected "
            << PART_ORD_UNIVERSAL << ", got: "
            << part_universal.mesh_meta_data_ordinal());

    ThrowRequireMsg(PART_ORD_OWNED == part_owned.mesh_meta_data_ordinal(),
            "Owned part ordinal is wrong, expected "
            << PART_ORD_OWNED << ", got: "
            << part_owned.mesh_meta_data_ordinal());

    ThrowRequireMsg(PART_ORD_SHARED == part_shared.mesh_meta_data_ordinal(),
            "Shared part ordinal is wrong, expected "
            << PART_ORD_SHARED << ", got: "
            << part_shared.mesh_meta_data_ordinal());

    //  Shared entities may have been modified due to relationship changes.
    //  Send just the current induced memberships from the sharing to
    //  the owning processes.
    {
        CommAll comm(p_comm);

        pack_induced_memberships(*this, comm, m_entity_comm_list);

        BABBLE_STK_PARALLEL_COMM(p_comm, "          internal_resolve_shared_membership calling allocate_buffers");
        comm.allocate_buffers(p_size / 4);

        pack_induced_memberships(*this, comm, m_entity_comm_list);

        BABBLE_STK_PARALLEL_COMM(p_comm, "          internal_resolve_shared_membership calling communicate");
        comm.communicate();

        OrdinalVector empty, induced_parts, current_parts, remove_parts;
        PartVector inducedParts, removeParts;

        for(EntityCommListInfoVector::iterator i = m_entity_comm_list.begin(); i != m_entity_comm_list.end(); ++i)
        {
            bool i_own_this_entity_in_comm_list = i->owner == p_rank;
            if( i_own_this_entity_in_comm_list )
            {
                // Receiving from all sharing processes

                empty.clear();
                induced_parts.clear();
                current_parts.clear();
                remove_parts.clear();

                induced_part_membership(*this, i->entity, empty, induced_parts);

                for(PairIterEntityComm ec = entity_comm_map_shared(i->key); !ec.empty(); ++ec)
                {
                    CommBuffer & buf = comm.recv_buffer(ec->proc);

                    unsigned count = 0;
                    buf.unpack<unsigned>(count);
                    for(unsigned j = 0; j < count; ++j)
                    {
                        unsigned part_ord = 0;
                        buf.unpack<unsigned>(part_ord);
                        insert_ordinal(induced_parts, part_ord);
                    }
                }

                // Remove any part that is an induced part but is not
                // in the induced parts list.

                this->bucket(i->entity).supersets(current_parts);

                OrdinalVector::const_iterator induced_parts_begin = induced_parts.begin(),
                        induced_parts_end = induced_parts.end();

                for(OrdinalVector::iterator
                p = current_parts.begin(); p != current_parts.end(); ++p)
                {
                    if(meta.get_parts()[*p]->was_induced(i->key.rank()) &&
                            !contains_ordinal(induced_parts_begin, induced_parts_end, *p))
                    {
                        remove_parts.push_back(*p);
                    }
                }

                inducedParts.clear();
                removeParts.clear();

                inducedParts.reserve(induced_parts.size());
                for(unsigned ipart = 0; ipart < induced_parts.size(); ++ipart)
                {
                    inducedParts.push_back(&m_mesh_meta_data.get_part(induced_parts[ipart]));
                }
                removeParts.reserve(remove_parts.size());
                for(unsigned ipart = 0; ipart < remove_parts.size(); ++ipart)
                {
                    removeParts.push_back(&m_mesh_meta_data.get_part(remove_parts[ipart]));
                }

                internal_change_entity_parts(i->entity, inducedParts, removeParts);
            }
        }
    }

    //------------------------------
    // The owners have complete knowledge of memberships.
    // Send membership information to sync the shared and ghosted copies.
    // Only need to do this for entities that have actually changed.

    {
        std::vector<EntityProc> send_list;

        generate_send_list(*this, m_sync_count, p_rank, m_entity_comm_list, send_list);

        CommAll comm(p_comm);

        pack_part_memberships(*this, comm, send_list);

        BABBLE_STK_PARALLEL_COMM(p_comm, "          internal_resolve_shared_membership calling allocate_buffers");
        comm.allocate_buffers(p_size / 4);

        pack_part_memberships(*this, comm, send_list);

        BABBLE_STK_PARALLEL_COMM(p_comm, "          internal_resolve_shared_membership calling communicate");
        comm.communicate();

        for(int p = 0; p < p_size; ++p)
        {
            CommBuffer & buf = comm.recv_buffer(p);
            while(buf.remaining())
            {

                PartVector owner_parts, current_parts, remove_parts;

                EntityKey key;
                buf.unpack<EntityKey>(key);
                unsigned count = 0;
                buf.unpack<unsigned>(count);
                for(unsigned j = 0; j < count; ++j)
                {
                    unsigned part_ord = 0;
                    buf.unpack<unsigned>(part_ord);
                    insert(owner_parts, *all_parts[part_ord]);
                }

                // Any current part that is not a member of owners_parts
                // must be removed.

                Entity const entity = find_entity(*this, m_entity_comm_list, key).entity;

                this->bucket(entity).supersets(current_parts);

                for(PartVector::iterator
                ip = current_parts.begin(); ip != current_parts.end(); ++ip)
                {
                    Part * const part = *ip;
                    const unsigned part_ord = part->mesh_meta_data_ordinal();
                    if(PART_ORD_UNIVERSAL != part_ord &&
                            PART_ORD_OWNED != part_ord &&
                            PART_ORD_SHARED != part_ord &&
                            !contain(m_ghost_parts, *part) &&
                            !contain(owner_parts, *part))
                    {
                        remove_parts.push_back(part);
                    }
                }

                internal_change_entity_parts(entity, owner_parts, remove_parts);
            }
        }
    }
}

void BulkData::internal_update_fast_comm_maps()
{
  if (parallel_size() > 1) {
    EntityCommListInfoVector const& all_comm = comm_list();

    // Flush previous map
    const EntityRank num_ranks = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());
    m_volatile_fast_shared_comm_map.resize(num_ranks);
    for (EntityRank r = stk::topology::BEGIN_RANK; r < num_ranks; ++r) {
      m_volatile_fast_shared_comm_map[r].resize(parallel_size());
      for (int proc = 0; proc < parallel_size(); ++proc) {
        m_volatile_fast_shared_comm_map[r][proc].clear();
      }
    }

    // Assemble map, find all shared entities and pack into volatile fast map
    for (size_t i = 0, ie = all_comm.size(); i < ie; ++i) {
      Entity const e        = all_comm[i].entity;
      EntityKey const key   = all_comm[i].key;
      MeshIndex const& idx  = mesh_index(e);
      EntityRank const rank = key.rank();

      FastMeshIndex fast_idx;
      fast_idx.bucket_id  = idx.bucket->bucket_id();
      fast_idx.bucket_ord = idx.bucket_ordinal;

      PairIterEntityComm ec = entity_comm_map(key);
      for (; !ec.empty() && ec->ghost_id == 0; ++ec) {
        m_volatile_fast_shared_comm_map[rank][ec->proc].push_back(fast_idx);
      }
    }

    // Need to shrink-to-fit these vectors?
  }
}



namespace impl {

unsigned get_ordinal(const Part* part)
{ return part->mesh_meta_data_ordinal(); }

unsigned get_ordinal(unsigned ord)
{ return ord; }

const Part& get_part(const Part* part, MetaData& meta)
{ return *part; }

const Part& get_part(unsigned ord, MetaData& meta)
{ return *meta.get_parts()[ord]; }

void filter_out( std::vector<unsigned> & vec ,
                 const std::vector<Part*> & parts ,
                 std::vector<Part*> & removed )
{
  std::vector<unsigned>::iterator i , j ;
  i = j = vec.begin();

  std::vector<Part*>::const_iterator ip = parts.begin() ;

  while ( j != vec.end() && ip != parts.end() ) {
    if      ( get_ordinal(*ip) < *j ) { ++ip ; }
    else if ( *j < get_ordinal(*ip) ) { *i = *j ; ++i ; ++j ; }
    else {
      removed.push_back( *ip );
      ++j ;
      ++ip ;
    }
  }

  if ( i != j ) { vec.erase( i , j ); }
}

void merge_in( std::vector<unsigned> & vec , const std::vector<Part*> & parts )
{
  std::vector<unsigned>::iterator i = vec.begin();
  std::vector<Part*>::const_iterator ip = parts.begin() ;

  for ( ; i != vec.end() && ip != parts.end() ; ++i ) {

    const unsigned ord = get_ordinal(*ip);

    if ( ord <= *i ) {
      if ( ord < *i ) { i = vec.insert( i , ord ); }
      // Now have: ord == *i
      ++ip ;
    }
  }

  for ( ; ip != parts.end() ; ++ip ) {
    const unsigned ord = get_ordinal(*ip);
    vec.push_back( ord );
  }
}

} // namespace impl

void BulkData::change_entity_parts( Entity entity,
    const PartVector & add_parts ,
    const PartVector & remove_parts,
    bool always_propagate_internal_changes)
{
    INCREMENT_ENTITY_MODIFICATION_COUNTER(PUBLIC, entity_rank(entity), DESTROY_ALL_GHOSTING);
    internal_verify_and_change_entity_parts(entity, add_parts, remove_parts, always_propagate_internal_changes);
}

void BulkData::internal_verify_and_change_entity_parts( Entity entity,
                                                        const PartVector & add_parts ,
                                                        const PartVector & remove_parts,
                                                        bool always_propagate_internal_changes)
{
  TraceIfWatching("stk::mesh::BulkData::internal_change_entity_parts", LOG_ENTITY, entity_key(entity));
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "entity state: " << entity_key(entity));

  require_ok_to_modify();

// When stk parallel is used within Fmwk, this assertion (require_entity_owner) is violated
// So, temporarily, don't test this assertion if SIERRA_MIGRATION is defined, and the bulk
// data point is set.  (Any other use case will go ahead and test this assertion.)
#ifdef SIERRA_MIGRATION
  if (NULL == get_fmwk_bulk_data())
    require_entity_owner( entity , m_parallel_rank );
#else
  require_entity_owner( entity , m_parallel_rank );
#endif

  const EntityRank ent_rank = entity_rank(entity);

  const EntityRank undef_rank  = InvalidEntityRank;

  // Transitive addition and removal:
  // 1) Include supersets of add_parts
  // 2) Do not include a remove_part if it appears in the add_parts
  // 3) Include subsets of remove_parts

  // most parts will at least have universal and topology part as supersets
  const unsigned expected_min_num_supersets = 2;

  PartVector a_parts;
  a_parts.reserve( std::distance(add_parts.begin(), add_parts.end()) * (expected_min_num_supersets + 1) );
  for(PartVector::const_iterator add_iter=add_parts.begin(); add_iter!=add_parts.end(); ++add_iter) {
#ifdef FMWK_NO_GLOBALLY_SHARED_ELEMENTS
    ThrowErrorMsgIf(entity_rank == stk::topology::ELEMENT_RANK && **add_iter == mesh_meta_data().globally_shared_part(), "FMWK_NO_GLOBALLY_SHARED_ELEMENTS  Error in BulkData::internal_change_entity_parts, trying to make an element globally shared!");
#endif // FMWK_NO_GLOBALLY_SHARED_ELEMENTS
    a_parts.push_back((*add_iter));
  }
  bool quick_verify_check = true;

  for ( PartVector::const_iterator ia = add_parts.begin(); ia != add_parts.end() ; ++ia ) {
    quick_verify_check = quick_verify_check &&
      internal_quick_verify_change_part(*ia, ent_rank, undef_rank);
    const PartVector& supersets = (*ia)->supersets();
    for(PartVector::const_iterator s_iter=supersets.begin(), s_end=supersets.end();
        s_iter!=s_end; ++s_iter) {
      a_parts.push_back((*s_iter));
    }
  }

  order(a_parts);

  PartVector::const_iterator a_parts_begin = a_parts.begin(),
                                a_parts_end   = a_parts.end();
  PartVector r_parts ;

  for ( PartVector::const_iterator ir = remove_parts.begin(); ir != remove_parts.end() ; ++ir ) {

    // The following guards should be in the public interface to
    // changing parts.  However, internal mechanisms such as changing
    // ownership calls this function to add or remove an entity from
    // the three special parts.  Without refactoring, these guards
    // cannot be put in place.
    /*
    ThrowErrorMsgIf( m_mesh_meta_data.universal_part() == **ir,
                     "Cannot remove entity from universal part" );
    ThrowErrorMsgIf( m_mesh_meta_data.locally_owned_part() == **ir,
                     "Cannot remove entity from locally owned part" );
    ThrowErrorMsgIf( m_mesh_meta_data.globally_shared_part() == **ir,
                     "Cannot remove entity from globally shared part" );
    */

    quick_verify_check = quick_verify_check &&
      internal_quick_verify_change_part(*ir, ent_rank, undef_rank);

    if ( ! contains_ordinal_part( a_parts_begin, a_parts_end , (*ir)->mesh_meta_data_ordinal() ) ) {
      r_parts.push_back( (*ir) );
      for ( PartVector::const_iterator  cur_part = (*ir)->subsets().begin() ;
            cur_part != (*ir)->subsets().end() ;
            ++cur_part )
        if ( bucket(entity).member ( **cur_part ) )
          r_parts.push_back ( (*cur_part) );
    }
  }

  order(r_parts);

  // If it looks like we have a problem, run the full check and we should
  // expect to see an exception thrown; otherwise, only do the full check in
  // debug mode because it incurs significant overhead.
  if ( ! quick_verify_check ) {
    internal_verify_change_parts( m_mesh_meta_data , entity , a_parts );
    internal_verify_change_parts( m_mesh_meta_data , entity , r_parts );
    ThrowRequireMsg(false, "Expected throw from verify methods above.");
  }
  else {
#ifndef NDEBUG
    internal_verify_change_parts( m_mesh_meta_data , entity , a_parts );
    internal_verify_change_parts( m_mesh_meta_data , entity , r_parts );
#endif
  }

  internal_change_entity_parts( entity , a_parts , r_parts , always_propagate_internal_changes );
}

//  The 'add_parts' and 'remove_parts' are complete and disjoint.
//  Changes need to have parallel resolution during
//  modification_end.

void BulkData::internal_change_entity_parts(
  Entity entity ,
  const std::vector<Part*> & add_parts ,
  const std::vector<Part*> & remove_parts,
  bool always_propagate_internal_changes )
{
  TraceIfWatching("stk::mesh::BulkData::internal_change_entity_parts", LOG_ENTITY, entity_key(entity));
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "entity state: " << entity_key(entity));
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "add_parts: " << add_parts);
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "remove_parts: " << remove_parts);

  Bucket * const bucket_old = bucket_ptr( entity );

  if ( bucket_old && bucket_old->member_all( add_parts ) &&
              ! bucket_old->member_any( remove_parts ) ) {
    // Is already a member of all add_parts,
    // is not a member of any remove_parts,
    // thus nothing to do.
    return ;
  }

  INCREMENT_ENTITY_MODIFICATION_COUNTER(INTERNAL, entity_rank(entity), CHANGE_ENTITY_PARTS);

  const unsigned locally_owned_ordinal = m_mesh_meta_data.locally_owned_part().mesh_meta_data_ordinal();

  bool add_to_locally_owned = false;
  for (std::vector<Part*>::const_iterator itr = add_parts.begin(), end_itr = add_parts.end(); itr != end_itr; ++itr) {
    if ( impl::get_ordinal(*itr) == locally_owned_ordinal ) {
      add_to_locally_owned = true;
      break;
    }
  }
  add_to_locally_owned = add_to_locally_owned && (!bucket_old || !bucket_old->owned());


  bool remove_from_locally_owned = false;
  for (std::vector<Part*>::const_iterator itr = remove_parts.begin(), end_itr = remove_parts.end(); itr != end_itr; ++itr) {
    if ( impl::get_ordinal(*itr) == locally_owned_ordinal ) {
      remove_from_locally_owned = true;
      break;
    }
  }
  remove_from_locally_owned = remove_from_locally_owned && (!bucket_old || bucket_old->owned());

  if (add_to_locally_owned) {

    ++m_closure_count[entity.local_offset()];

    // update downward connectivity closure count
    if (bucket_old) {
      for (EntityRank rank = stk::topology::NODE_RANK, end_rank = bucket_old->entity_rank(); rank < end_rank; ++rank) {
        unsigned num = num_connectivity(entity,rank);
        Entity const * entities = begin(entity,rank);
        for (unsigned i =0; i<num; ++i) {
          ++m_closure_count[entities[i].local_offset()];
        }
      }
    }

  }
  else if (remove_from_locally_owned)
  {
    --m_closure_count[entity.local_offset()];

    // update downward connectivity closure count
    if (bucket_old) {
      for (EntityRank rank = stk::topology::NODE_RANK, end_rank = bucket_old->entity_rank(); rank < end_rank; ++rank) {
        unsigned num = num_connectivity(entity,rank);
        Entity const * entities = begin(entity,rank);
        for (unsigned i =0; i<num; ++i) {
          --m_closure_count[entities[i].local_offset()];
        }
      }
    }
  }

  std::vector<Part*> parts_removed ;

  OrdinalVector parts_total ; // The final part list

  //--------------------------------



  if ( bucket_old ) {
    // Keep any of the existing bucket's parts
    // that are not a remove part.
    // This will include the 'intersection' parts.
    //
    // These parts are properly ordered and unique.

    const std::pair<const unsigned *, const unsigned*>
      bucket_parts = bucket_old->superset_part_ordinals();

    const unsigned * parts_begin = bucket_parts.first;
    const unsigned * parts_end   = bucket_parts.second;

    const unsigned num_bucket_parts = parts_end - parts_begin;
    parts_total.reserve( num_bucket_parts + add_parts.size() );
    parts_total.insert( parts_total.begin(), parts_begin , parts_end);

    if ( !remove_parts.empty() ) {
      parts_removed.reserve(remove_parts.size());
      impl::filter_out( parts_total , remove_parts , parts_removed );
    }
  }
  else {
    parts_total.reserve(add_parts.size());
  }

  if ( !add_parts.empty() ) {
    impl::merge_in( parts_total , add_parts );
  }

  if ( parts_total.empty() ) {
    // Always a member of the universal part.
    const unsigned univ_ord =
      m_mesh_meta_data.universal_part().mesh_meta_data_ordinal();
    parts_total.push_back( univ_ord );
  }

  EntityRank e_rank = entity_rank(entity);

  //--------------------------------
  // Move the entity to the new partition.
  stk::mesh::impl::Partition *partition =
          m_bucket_repository.get_or_create_partition(e_rank, parts_total);

  if (bucket_old) {
    bucket_old->getPartition()->move_to(entity, *partition);
  }
  else {
    partition->add(entity);
  }

  ////
  //// SHOULD WE FIND A WAY TO PUT THE REST OF THIS IN Partition::move_to(..)?
  ////

  // Propagate part changes through the entity's relations.
  //(Only propagate part changes for parts which have a primary-entity-rank that matches
  // the entity's rank. Other parts don't get induced...)

  std::vector<Part*> inducable_parts_removed;
  for(std::vector<Part*>::const_iterator pr=parts_removed.begin(), prend=parts_removed.end(); pr!=prend; ++pr) {
    Part const& check_part = impl::get_part(*pr, m_mesh_meta_data);
    if (check_part.should_induce(e_rank)) {
      inducable_parts_removed.push_back(*pr);
    }
  }

  if (always_propagate_internal_changes ||
      !inducable_parts_removed.empty() ) {
    internal_propagate_part_changes( entity , inducable_parts_removed );
  }
}


//----------------------------------------------------------------------
// Deduce propagation of part membership changes to a 'from' entity
// to the related 'to' entities.  There can be both additions and
// removals.

void BulkData::internal_propagate_part_changes(
  Entity entity ,
  const std::vector<Part*> & removed )
{
  TraceIfWatching("stk::mesh::BulkData::internal_propagate_part_changes",
      LOG_ENTITY,
      entity_key(entity));
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "entity state: " << entity_key(entity));
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "Removed: " << removed);

  m_check_invalid_rels = false;

  const EntityRank erank = entity_rank(entity);
  const EntityRank end_rank = static_cast<EntityRank>(m_mesh_meta_data.entity_rank_count());

  OrdinalVector to_del , to_add , empty ;
  EntityVector temp_entities;
  for (EntityRank irank = stk::topology::BEGIN_RANK; irank < erank; ++irank)
  {
    size_t num_rels = num_connectivity(entity, irank);
    if (num_rels > 0) {
      Entity const *rel_entities = begin(entity, irank);
      for (size_t j = 0; j < num_rels; ++j)
      {
        Entity e_to = rel_entities[j];

        if (e_to == Entity::InvalidEntity)
        {
          continue;
        }

        to_del.clear();
        to_add.clear();
        empty.clear();

        // Induce part membership from this relationship to
        // pick up any additions.
        induced_part_membership(*this, entity, empty, irank, to_add );

        if ( ! removed.empty() ) {
          // Something was removed from the 'from' entity,
          // deduce what may have to be removed from the 'to' entity.

          // Deduce parts for 'e_to' from all upward relations.
          // Any non-parallel part that I removed that is not deduced for
          // 'e_to' must be removed from 'e_to'

          EntityRank e_to_rank = entity_rank(e_to);

          Entity const* back_rel_entities = NULL;
          int num_back_rels = 0;
          for (EntityRank to_rel_rank_i = static_cast<EntityRank>(e_to_rank + 1); to_rel_rank_i < end_rank; ++to_rel_rank_i)
          {
            if (connectivity_map().valid(e_to_rank, to_rel_rank_i)) {
              num_back_rels     = num_connectivity(e_to, to_rel_rank_i);
              back_rel_entities = begin(e_to, to_rel_rank_i);
            }
            else {
              num_back_rels = get_connectivity(*this, e_to, to_rel_rank_i, temp_entities);
              back_rel_entities = &*temp_entities.begin();
            }

            for (int k = 0; k < num_back_rels; ++k)
            {
              if (entity != back_rel_entities[k])  // Already did this entity
              {
                // Relation from to_rel->entity() to e_to
                induced_part_membership(*this, back_rel_entities[k], empty, e_to_rank, to_add );
              }
            }
          }

          OrdinalVector::const_iterator to_add_begin = to_add.begin(),
            to_add_end   = to_add.end();

          for (std::vector<Part*>::const_iterator
                  k = removed.begin() ; k != removed.end() ; ++k ) {
            if ( ! contains_ordinal( to_add_begin, to_add_end , impl::get_ordinal(*k) ) ) {
              induced_part_membership( impl::get_part(*k, m_mesh_meta_data), erank, irank, to_del );
            }
          }
        }


        PartVector addParts, delParts, emptyParts;

        addParts.reserve(to_add.size());
        for(unsigned ipart=0; ipart<to_add.size(); ++ipart) {
          addParts.push_back(&m_mesh_meta_data.get_part(to_add[ipart]));
        }
        delParts.reserve(to_del.size());
        for(unsigned ipart=0; ipart<to_del.size(); ++ipart) {
          delParts.push_back(&m_mesh_meta_data.get_part(to_del[ipart]));
        }



        if ( parallel_size() < 2 || entity_comm_map_shared(entity_key(e_to)).empty() ) {
          // Entirely local, ok to remove memberships now
          internal_change_entity_parts( e_to , addParts , delParts );
        }
        else {
          // Shared, do not remove memberships now.
          // Wait until modification_end.
          internal_change_entity_parts( e_to , addParts , emptyParts );
        }
      }
    }
  }
  m_check_invalid_rels = true;
}

// TODO Change the methods below to requirements (private, const invariant checkers)

// Do not allow any of the induced part memberships to explicitly
// appear in the add or remove parts lists.
// 1) Intersection part
// 3) Part that does not match the entity rank.

void BulkData::internal_verify_change_parts( const MetaData   & meta ,
                                             const Entity entity ,
                                             const std::vector<Part*> & parts ) const
{
  const std::vector<std::string> & rank_names = meta.entity_rank_names();
  const EntityRank undef_rank  = InvalidEntityRank;
  const EntityRank erank = entity_rank(entity);

  bool ok = true ;
  std::ostringstream msg ;

  for (std::vector<Part*>::const_iterator
        i = parts.begin() ; i != parts.end() ; ++i ) {

    const Part & p = impl::get_part(*i, m_mesh_meta_data);
    const unsigned part_rank = p.primary_entity_rank();

    bool intersection_ok=false, rel_target_ok=false, rank_ok=false;
    internal_basic_part_check(&p, erank, undef_rank, intersection_ok, rel_target_ok, rank_ok);

    if ( !intersection_ok || !rel_target_ok || !rank_ok ) {
      if ( ok ) {
        ok = false ;
        msg << "change parts for entity " << identifier(entity);
        msg << " , { " ;
      }
      else {
        msg << " , " ;
      }

      msg << p.name() << "[" ;
      if ( part_rank < rank_names.size() ) {
        msg << rank_names[ part_rank ];
      }
      else {
        msg << part_rank ;
      }
      msg << "] " ;
      if ( !intersection_ok ) { msg << "is_intersection " ; }
      if ( !rel_target_ok )   { msg << "is_relation_target " ; }
      if ( !rank_ok )         { msg << "is_bad_rank " ; }
    }
  }

  ThrowErrorMsgIf( !ok, msg.str() << "}" );
}

} // namespace mesh
} // namespace stk
