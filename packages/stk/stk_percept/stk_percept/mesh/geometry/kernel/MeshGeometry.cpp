#include "MeshGeometry.hpp"

MeshGeometry::MeshGeometry(GeometryKernel* geom)
{
  geomKernel = geom;
  mDbgNodeCoords[0] = -0.00477133907617983;
  mDbgNodeCoords[1] = -0.00477133907617983;
  mDbgNodeCoords[2] =  0.260484055257467;
}

MeshGeometry::~MeshGeometry()
{

}

void MeshGeometry::add_evaluator(GeometryEvaluator* evaluator)
{
  geomEvaluators.push_back(evaluator);
}

void MeshGeometry::add_evaluators(std::vector<GeometryEvaluator*> evaluators)
{
  geomEvaluators.insert(geomEvaluators.end(), evaluators.begin(), evaluators.end());
  // std::for_each( evaluators.begin(),
  //                evaluators.end(),
  //                std::mem_fun( &MeshGeometry::add_evaluator ) );
}

const std::vector<GeometryEvaluator*>& MeshGeometry::getGeomEvaluators()
{
  return geomEvaluators;
}

/**
 * Return 0,1,2,3 if the node is on a geometry vertex, curve, surface or domain.
 */
int MeshGeometry::classify_node(const stk::mesh::Entity& node, std::vector<size_t>& curveEvaluators, std::vector<size_t>& surfEvaluators)
{
  const stk::mesh::Bucket& bucket = node.bucket();
  return classify_bucket(bucket, curveEvaluators, surfEvaluators);
}

/**
 * Return 0,1,2,3 if the bucket is on a geometry vertex, curve, surface or domain.
 */
int MeshGeometry::classify_bucket(const stk::mesh::Bucket& bucket, std::vector<size_t>& curveEvaluators, std::vector<size_t>& surfEvaluators)
{
  // Each bucket contains the set of nodes with unique part intersections.
  // This means that every nodes will be in exactly one bucket.  But, the
  // nodes on curves are also in the part for the adjacent surfaces which
  // means more than one evaluator will be selected for those buckets which
  // are on the boundary of an entity (i.e. buckets representing curves
  // and vertices).  We first create a list of all the evaluators that
  // might be relevant.

  curveEvaluators.resize(0);
  surfEvaluators.resize(0);

  size_t s;
  for (s=0; s<geomEvaluators.size(); s++)
    {
      if (geomEvaluators[s]->mMesh(bucket))
        {
          if (geomKernel->is_curve(s))
            {
              curveEvaluators.push_back(s);
            }
          else if (geomKernel->is_surface(s))
            {
              surfEvaluators.push_back(s);
            }
        }
    }

  if ( curveEvaluators.size() > 1 )
    {
      // This is a bucket representing a vertex.  No need to do anything
      // since the node will already be on the vertex, and no new nodes
      // ever get created assigned to vertices during refinement.
      //std::cout << "Vertex node encountered" << std::endl;
      return 0;
    }
  if ( curveEvaluators.size() == 1 )
    {
      // This bucket represents a geometric curve.  Snap to it.
      //std::cout << "Snapping to curve" << curveEvaluators[0] << std::endl;
      //snap_nodes( mesh_data, bucket, curveEvaluators[0] );
      return 1;
    }

  if ( surfEvaluators.size() == 1 )
    {
      //std::cout << "Snapping to surface" << surfEvaluators[0] << std::endl;
      // This bucket represents a geometric surface.  Snap to it.
      //snap_nodes( mesh_data, bucket, surfEvaluators[0] );
      return 2;
    }
  return -1;  // error condition
}

void MeshGeometry::snap_points_to_geometry(PerceptMesh* mesh_data)
{
  BulkData& bulkData = *mesh_data->getBulkData();

  const std::vector<Bucket*> & buckets = bulkData.buckets( stk::mesh::fem::FEMMetaData::NODE_RANK );

  for ( std::vector<Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
  {
    Bucket& bucket = **k;
    if ( contains_dbg_node( mesh_data, bucket ) )
    {
      std::cout << "     DBG Node FOUND" << std::endl;
    }

    std::vector<size_t> curveEvaluators;
    std::vector<size_t> surfEvaluators;

    int type = classify_bucket(bucket, curveEvaluators, surfEvaluators);
    switch (type) {
    case 0:
      // This is a bucket representing a vertex.  No need to do anything
      // since the node will already be on the vertex, and no new nodes
      // ever get created assigned to vertices during refinement.
      //std::cout << "Vertex node encountered" << std::endl;
      break;
    case 1:
      // This bucket represents a geometric curve.  Snap to it.
      //std::cout << "Snapping to curve" << curveEvaluators[0] << std::endl;
      snap_nodes( mesh_data, bucket, curveEvaluators[0] );
      break;
    case 2:
      //std::cout << "Snapping to surface" << surfEvaluators[0] << std::endl;
      // This bucket represents a geometric surface.  Snap to it.
      snap_nodes( mesh_data, bucket, surfEvaluators[0] );
      break;
    case -1:
    default:
      printf( "ERROR: A bucket found without a geometric evaluator.\n" );
      break;
    }


  }
}

void MeshGeometry::snap_points_to_geometry(PerceptMesh* mesh_data, std::vector<stk::mesh::Entity *>& nodes)
{
  for (unsigned inode=0; inode < nodes.size(); inode++)
  {
    Bucket& bucket = nodes[inode]->bucket();

    if ( contains_dbg_node( mesh_data, bucket ) )
    {
      std::cout << "     DBG Node FOUND" << std::endl;
    }

    // Each bucket contains the set of nodes with unique part intersections.
    // This means that every nodes will be in exactly one bucket.  But, the
    // nodes on curves are also in the part for the adjacent surfaces which
    // means more than one evaluator will be selected for those buckets which
    // are on the boundary of an entity (i.e. buckets representing curves
    // and vertices).  We first create a list of all the evaluators that
    // might be relevant.
    std::vector<size_t> curveEvaluators;
    std::vector<size_t> surfEvaluators;

    int type = classify_bucket(bucket, curveEvaluators, surfEvaluators);
    switch (type) {
    case 0:
      // This is a bucket representing a vertex.  No need to do anything
      // since the node will already be on the vertex, and no new nodes
      // ever get created assigned to vertices during refinement.
      //std::cout << "Vertex node encountered" << std::endl;
      break;
    case 1:
      // This bucket represents a geometric curve.  Snap to it.
      //std::cout << "Snapping to curve" << curveEvaluators[0] << std::endl;
      snap_node( mesh_data, *nodes[inode], curveEvaluators[0] );
      break;
    case 2:
      //std::cout << "Snapping to surface" << surfEvaluators[0] << std::endl;
      // This bucket represents a geometric surface.  Snap to it.
      snap_node( mesh_data, *nodes[inode], surfEvaluators[0] );
      break;
    case -1:
    default:
      printf( "ERROR: A bucket found without a geometric evaluator.\n" );
      break;
    }


  }
}

void MeshGeometry::snap_node
(
  PerceptMesh *mesh_data,
  Entity & node,
  size_t evaluator_idx
)
{
  VectorFieldType* coordField = mesh_data->getCoordinatesField();

  Part* new_nodes_part = mesh_data->getNonConstPart("refine_new_nodes_part");
  Selector new_nodes_part_selector;
  if (new_nodes_part) new_nodes_part_selector = Selector(*new_nodes_part);

  std::string str = geomKernel->get_attribute(evaluator_idx);
  {

    double * coord = stk::mesh::field_data( *coordField , node );
    double delta[3] = {coord[0], coord[1], coord[2]};
    bool doPrint = DEBUG_GEOM_SNAP && new_nodes_part_selector(node);
    if (doPrint)
    {
      std::cout << "tmp geom snap_points_to_geometry eval name= " << str << " node id= " << node.identifier() 
                << " coords b4= " << coord[0] << " " << coord[1] << " " << coord[2];
    }

    if ( is_dbg_node( coord ) )
    {
      std::cout << "Node in question being projected" << std::endl;
    }

    geomKernel->snap_to(coord, geomEvaluators[evaluator_idx]->mGeometry);

    if (doPrint)
    {
      delta[0] = coord[0] - delta[0];
      delta[1] = coord[1] - delta[1];
      delta[2] = coord[2] - delta[2];
      double dtot = std::sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]);
      std::cout << " coords af= " << coord[0] << " " << coord[1] << " " << coord[2] 
                << " delta= " << delta[0] << " " << delta[1] << " " << delta[2] 
                << " deltaTot= " << dtot
                << " " << str.substr(6,5)
                << std::endl;
      //if (str.substr(6,5)=="20004") block_20004
      //{
      //  std::cout << "found 20004" << std::endl;
      //}
    }
  }
}

void MeshGeometry::snap_nodes
(
  PerceptMesh *mesh_data,
  Bucket &bucket,
  size_t evaluator_idx
)
{
  //VectorFieldType* coordField = mesh_data->getCoordinatesField();
  const unsigned num_nodes_in_bucket = bucket.size();

  Part* new_nodes_part = mesh_data->getNonConstPart("refine_new_nodes_part");
  Selector new_nodes_part_selector;
  if (new_nodes_part) new_nodes_part_selector = Selector(*new_nodes_part);

  std::string str = geomKernel->get_attribute(evaluator_idx);
  for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
  {
    Entity& node = bucket[iNode];

    snap_node(mesh_data, node, evaluator_idx);
  }
}

bool MeshGeometry::contains_dbg_node
(
  PerceptMesh *mesh_data,
  Bucket &bucket
)
{
  VectorFieldType* coordField = mesh_data->getCoordinatesField();
  const unsigned num_nodes_in_bucket = bucket.size();

  for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
  {
    Entity& node = bucket[iNode];

    double * coord = stk::mesh::field_data( *coordField , node );
    if ( is_dbg_node( coord ) )
    {
      return true;
    }
  }
  return false;
}

bool MeshGeometry::is_dbg_node( double node_coord[3] )
{
  double dx = node_coord[0] - mDbgNodeCoords[0];
  double dy = node_coord[1] - mDbgNodeCoords[1];
  double dz = node_coord[2] - mDbgNodeCoords[2];
  double dist = sqrt( (dx*dx) + (dy*dy) + (dz*dz) );
  if ( dist < 0.0001 )
  {
    return true;
  }
  return false;
}

