#if 0
/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <vector>
#include <fstream>

#include <assert.h>

#include <stk_io/StkMeshIoBroker.hpp>
#include <init/Ionit_Initializer.h>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/WriterExt.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelIndex.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search_util/stk_mesh/PrintEntityProc.hpp>
#include <transfer/UseCaseIsInElement.hpp>

#include <stk_search_util/stk_mesh/CreateBoundingBox.hpp>
#include <stk_search_util/stk_mesh/PrintBoundingBox.hpp>

using namespace stk::diag;
using namespace use_case;

struct EntityKeyDecomp {
  unsigned operator()( const unsigned p_size ,
                       const stk::mesh::EntityKey & key ) const
    { return ( key.id() >> 8 ) % p_size ; }
};

typedef stk::util::ParallelIndex<stk::mesh::EntityKey,unsigned,EntityKeyDecomp> ParallelIndex;


namespace {

template <class T> static std::string to_string(const T & t)
{
  std::ostringstream os;
  os << t;
  return os.str();
}


template<typename T>
std::ostream& operator<<(std::ostream& ostr, const std::vector<T> &v){
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(ostr, ", "));
  return ostr;
}


void get_idents(stk::mesh::BulkData &bulk_data,
                std::vector<ParallelIndex::Key> &ident_vector)
{
  BOOST_STATIC_ASSERT(sizeof(ParallelIndex::Key) >= sizeof(stk::mesh::EntityKey));
  ident_vector.clear();
  const stk::mesh::MetaData&   meta_data = stk::mesh::MetaData::get(bulk_data);

  std::vector<stk::mesh::Entity> entities;
  stk::mesh::Selector selector = meta_data.locally_owned_part();
  get_selected_entities(selector, bulk_data.buckets(stk::mesh::MetaData::NODE_RANK), entities);
  const size_t num_entities = entities.size();

  for (size_t i = 0; i < num_entities; ++i) {
//    const ParallelIndex::Key   p = entities[i]->key().value();
    const ParallelIndex::Key   p = bulk_data.entity_key(entities[i]);
    ident_vector.push_back(p);
  }
}

void check_query( const std::vector<ParallelIndex::Key> &recv_global_id_vector,
                  const std::vector<ParallelIndex::KeyProc>   &processor_numbers)
{
  std::vector< ParallelIndex::Key>::const_iterator
    r_i = recv_global_id_vector.begin(),
    r_e = recv_global_id_vector.end();

  std::vector<ParallelIndex::KeyProc>::const_iterator
    p_i = processor_numbers.begin(),
    p_e = processor_numbers.end();

  for (; r_e != r_i && p_e != p_i; ++r_i, ++p_i) {
    const stk::mesh::EntityKey ri    = *r_i;
    const stk::mesh::EntityKey pi    = p_i->first;
    const stk::mesh::EntityKey pi_p1 = (p_e == p_i + 1) ? stk::mesh::EntityKey(pi.rank(), pi.id() + 1) : (p_i+1)->first;
    if (pi == pi_p1) {
      std::cout << pi;

      // TODO: These cerr statements should be changed to ThrowErrorMsgIf
      std::cerr
        << " A send mesh node with global id "<<pi
        << " was found multiple times on different processors.  Since only active and locally\n"
        << " owned nodes are searched, a unique object should have been found.  This node was found\n"
        << " on processor number "<<p_i->second<<" and on processor "<<(p_i+1)->second<<".\n"
        << std::endl << StackTrace;
      std::exit(EXIT_FAILURE);
    } else if (ri < pi) {
      std::cerr
        << " An active and locally owned receiving mesh node with global id "<<ri
        << " was to be paired with a sending node with the same global id.\n"
        << " After an exaustive search of all the sending mesh parts on all processors,"
        << " no node was found with the correct global id.\n"
        << std::endl << StackTrace;
      std::exit(EXIT_FAILURE);
    } else if (pi < ri) {
      std::cerr
        << " An active and locally owned receiving mesh node with global id "<<ri
        << " was to be paired with a sending node with the same global id.\n"
        << " The parallel global search instead returned an object with global id "<<pi
        << " which makes no sense.\n"
        << std::endl << StackTrace;
      std::exit(EXIT_FAILURE);
    }
  }
  if (r_e != r_i) {
    const stk::mesh::EntityKey ri = *r_i;
    std::cerr
      << " An active and locally owned receiving mesh node with global id "<<ri
      << " was to be paired with a sending node with the same global id.\n"
      << " After an exaustive search of all the sending mesh parts on all processors,"
      << " no node was found with the correct global id.\n"
      << std::endl << StackTrace;
    std::exit(EXIT_FAILURE);
  }
  if (p_e != p_i) {
    const stk::mesh::EntityKey pi    = p_i->first;
    const stk::mesh::EntityKey pi_p1 = (p_e == p_i + 1) ? stk::mesh::EntityKey(pi.rank(), pi.id() + 1) : (p_i+1)->first;
    if (pi == pi_p1) {
      std::cerr
        << " A send mesh node with global id "<<pi
        << " was found multiple times on different processors.  Since only active and locally\n"
        << " owned nodes are searched, a unique object should have been found.  This node was found\n"
        << " on processor number "<<p_i->second<<" and on processor "<<(p_i+1)->second<<".\n"
        << std::endl << StackTrace;
      std::exit(EXIT_FAILURE);
    } else  {
      std::cerr
        << " The parallel global search returned an object with global id "<<pi
        << " which makes no sense since no object with that id was searched for.\n"
        << std::endl << StackTrace;
      std::exit(EXIT_FAILURE);
    }
  }
}

}

void
use_case_0_driver(stk::ParallelMachine  comm,
                  const std::string &working_directory,
                  const std::string &range_mesh_filename,
                  const std::string &range_mesh_type,
                  const std::string &range_entity,
                  const std::string &domain_mesh_filename,
                  const std::string &domain_mesh_type,
                  const std::string &domain_entity)
{
  stk::diag::WriterThrowSafe _write_throw_safe(dw());

  stk::CommAll comm_all( comm );

  dw().m(LOG_TRANSFER) << "Use case 2: Point (range) Point (domain) Copy Search" << stk::diag::dendl;
  dw().m(LOG_TRANSFER) << "Range  Entity Type = " << range_entity  << stk::diag::dendl;
  dw().m(LOG_TRANSFER) << "Domain Entity Type = " << domain_entity << stk::diag::dendl;

  stk::io::StkMeshIoBroker range_mesh_data(comm);
  std::string filename = working_directory + range_mesh_filename;
  range_mesh_data.open_mesh_database(filename, range_mesh_type);
  range_mesh_data.create_input_mesh();
  stk::mesh::MetaData &range_meta_data = range_mesh_data.meta_data();
  range_meta_data.commit();

  range_mesh_data.populate_bulk_data();
  stk::mesh::BulkData &range_bulk_data = range_mesh_data.bulk_data();

  stk::io::StkMeshIoBroker domain_mesh_data(comm);
  filename = working_directory + domain_mesh_filename;
  domain_mesh_data.open_mesh_database(filename, domain_mesh_type);
  domain_mesh_data.create_input_mesh();
  stk::mesh::MetaData &domain_meta_data = domain_mesh_data.meta_data();
  domain_meta_data.commit();

  domain_mesh_data.populate_bulk_data();
  stk::mesh::BulkData &domain_bulk_data = domain_mesh_data.bulk_data();

  // For this use case, the domain consists of an axis-aligned
  // bounding box for each 'domain_entity' in the mesh.  The range is a
  // also an axis-aligned  bounding box.

  std::vector<PointBoundingBox3D> range_vector;

  std::vector< ParallelIndex::Key > range_global_id_vector;
  std::vector< ParallelIndex::Key > domain_global_id_vector;
  get_idents(range_bulk_data,   range_global_id_vector);
  get_idents(domain_bulk_data,  domain_global_id_vector);

  dw().m(LOG_TRANSFER) << "range  " << range_global_id_vector  << stk::diag::dendl;
  dw().m(LOG_TRANSFER) << "domain " << domain_global_id_vector << stk::diag::dendl;

  ParallelIndex parallel_index(comm, range_global_id_vector);
  std::vector<ParallelIndex::KeyProc> processor_numbers;
  parallel_index.query( domain_global_id_vector, processor_numbers);

  check_query( domain_global_id_vector, processor_numbers);

  for (std::vector<ParallelIndex::KeyProc>::const_iterator i=processor_numbers.begin();
       i != processor_numbers.end(); ++i)
  {
//     stk::mesh::EntityKey entity_key;
//     entity_key.value(i->first);
    stk::mesh::EntityKey entity_key(i->first);

    const unsigned entity_rank = entity_key.rank();
    const stk::mesh::EntityId entity_id =  entity_key.id();
    const std::string & entity_rank_name = domain_meta_data.entity_rank_name( entity_rank );
    dw().m(LOG_TRANSFER)<<" contains "<<" "<<entity_rank_name<<"["<<entity_id<<"] Proc:"<<i->second<<stk::diag::dendl;
  }
}
#endif
