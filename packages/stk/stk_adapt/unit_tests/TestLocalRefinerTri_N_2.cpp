#include <stdlib.h>

#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <unit_tests/TestLocalRefinerTri_N_2.hpp>

//#define STK_PERCEPT_HAS_GEOMETRY
#undef STK_PERCEPT_HAS_GEOMETRY
#if defined( STK_PERCEPT_HAS_GEOMETRY )
#include <stk_adapt/geometry/GeometryKernelOpenNURBS.hpp>
#include <stk_adapt/geometry/MeshGeometry.hpp>
#include <stk_adapt/geometry/GeometryFactory.hpp>
#endif

// FIXME
// #include <stk_mesh/baseImpl/EntityImpl.hpp>
// #include <stk_mesh/base/Entity.hpp>
// FIXME

namespace stk {
  namespace adapt {


    // This is a very specialized test that is used in unit testing only (see unit_localRefiner/break_tri_to_tri_N_2 in UnitTestLocalRefiner.cpp)

    TestLocalRefinerTri_N_2::TestLocalRefinerTri_N_2(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk::mesh::FieldBase *proc_rank_field) :
      Refiner(eMesh, bp, proc_rank_field)
    {
    }


    void TestLocalRefinerTri_N_2::
    refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element, vector<NeededEntityType>& needed_entity_ranks)
    {
      //static int n_seq = 400;

      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);

      CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::mesh::MetaData::NODE_RANK);

      VectorFieldType* coordField = m_eMesh.get_coordinates_field();

      for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
        {
          unsigned numSubDimNeededEntities = 0;
          stk::mesh::EntityRank needed_entity_rank = needed_entity_ranks[ineed_ent].first;

          if (needed_entity_rank == m_eMesh.edge_rank())
            {
              numSubDimNeededEntities = cell_topo_data->edge_count;
            }
          else if (needed_entity_rank == m_eMesh.face_rank())
            {
              numSubDimNeededEntities = cell_topo_data->side_count;
            }
          else if (needed_entity_rank == stk::mesh::MetaData::ELEMENT_RANK)
            {
              numSubDimNeededEntities = 1;
            }

          // see how many edges are already marked
          int num_marked=0;
          if (needed_entity_rank == m_eMesh.edge_rank())
            {
              for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                {
                  bool is_empty = m_nodeRegistry->is_empty( element, needed_entity_rank, iSubDimOrd);
                  if (!is_empty) ++num_marked;
                }
            }

          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              /// note: at this level of granularity we can do single edge refinement, hanging nodes, etc.
              //SubDimCell_SDSEntityType subDimEntity;
              //getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
              //bool is_empty = m_nodeRegistry->is_empty( element, needed_entity_rank, iSubDimOrd);
              //if(1||!is_empty)

              if (needed_entity_rank == m_eMesh.edge_rank())
                {
                  stk::mesh::Entity node0 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
                  stk::mesh::Entity node1 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
                  double * const coord0 = eMesh.field_data( *coordField , node0 );
                  double * const coord1 = eMesh.field_data( *coordField , node1 );

                  // vertical line position
                  const double vx = 0.21;

                  // choose to refine or not
                  if ( std::fabs(coord0[0]-coord1[0]) > 1.e-3 &&
                       ( (coord0[0] < vx && vx < coord1[0]) || (coord1[0] < vx && vx < coord0[0]) )
                       )
                    {
                      (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, true);
                    }
                }

            } // iSubDimOrd
        } // ineed_ent
    }

    ElementUnrefineCollection TestLocalRefinerTri_N_2::buildTestUnrefineList()
    {
      ElementUnrefineCollection elements_to_unref(*m_eMesh.get_bulk_data());

      const vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( stk::mesh::MetaData::ELEMENT_RANK );

      for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          //if (removePartSelector(**k))
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity element = bucket[ientity];

                // FIXME
                // skip elements that are already a parent (if there's no family tree yet, it's not a parent, so avoid throwing an error is isParentElement)
                const bool check_for_family_tree = false;
                bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);

                if (isParent)
                  continue;

                const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::mesh::MetaData::NODE_RANK);

                if (elem_nodes.size() && m_eMesh.isChildWithoutNieces(element, false) )
                  {
                    bool found = true;
                    for (unsigned inode=0; inode < elem_nodes.size(); inode++)
                      {
                        stk::mesh::Entity node = elem_nodes[inode].entity();
                        double *coord = eMesh.field_data( *m_eMesh.get_coordinates_field(), node );
                        if (coord[0] > 2.1 || coord[1] > 2.1)
                          {
                            found = false;
                            break;
                          }
                      }
                    if (found)
                      {
                        elements_to_unref.insert(element);
                        //std::cout << "tmp element id= " << m_eMesh.identifier(element) << " ";
                        //m_eMesh.print_entity(std::cout, element);
                      }
                  }
              }
          }
        }

      return elements_to_unref;
    }


  } // namespace adapt
} // namespace stk
