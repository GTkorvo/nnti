/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
/*!
  \file   PatchData.cpp

  \author Thomas Leurent
  \author Michael Brewer
  \date   2002-01-17
*/

#include "PatchData.hpp"
#include "MsqMeshEntity.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MeshInterface.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"
#include "GlobalPatch.hpp"
#include "PatchIterator.hpp"
#include "ExtraData.hpp"
#include "MappingFunctionSet.hpp"
#include "MappingFunction.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <list.h>
#  include <vector.h>
#  include <map.h>
#  include <algorithm.h>
#  include <numeric.h>
#  include <functional.h>
#  include <utility.h>
#else
#  include <list>
#  include <vector>
#  include <map>
#  include <algorithm>
#  include <numeric>
#  include <functional>
#  include <utility>
   using std::list;
   using std::map;
   using std::vector;
#endif

#ifdef MSQ_USE_OLD_IO_HEADERS
#  include <iostream.h>
#else
#  include <iostream>
   using std::ostream;
   using std::endl;
#endif

namespace Mesquite {

PatchData::PatchData()
  : myMesh(0),
    myDomain(0),
    numFreeVertices(0),
    numSlaveVertices(0),
    haveComputedInfos(0),
    dataList(0),
    mappingFuncs(0)
{
}


// Destructor
PatchData::~PatchData()
{
  notify_patch_destroyed();
}


void PatchData::get_minmax_element_unsigned_area(double& min, double& max, MsqError &err)
{
  if (!have_computed_info(MAX_UNSIGNED_AREA) ||
      !have_computed_info(MIN_UNSIGNED_AREA))
  {
    max=0;
    min=MSQ_DBL_MAX;
    size_t count = num_elements();
    for (size_t i=0; i<count; ++i) {
      double vol;
      vol = elementArray[i].compute_unsigned_area(*this, err); MSQ_ERRRTN(err);
      if (vol > max)
        max = vol;
      if (vol < min)
        min = vol;
    }
    note_have_info(MAX_UNSIGNED_AREA);
    note_have_info(MIN_UNSIGNED_AREA);
    computedInfos[MAX_UNSIGNED_AREA] = max;
    computedInfos[MIN_UNSIGNED_AREA] = min;
  }
  else
  {
    max = computedInfos[MAX_UNSIGNED_AREA];
    min = computedInfos[MIN_UNSIGNED_AREA];
  }
    
  if (max <= 0 || min < 0 || min == MSQ_DBL_MAX)
    MSQ_SETERR(err)(MsqError::INTERNAL_ERROR);
}

double PatchData::get_barrier_delta(MsqError &err)
{
  double result;
  if (have_computed_info(MINMAX_SIGNED_DET3D))
  {
    result = computedInfos[MINMAX_SIGNED_DET3D];
  }
  else
  {
    double min= MSQ_DBL_MAX;
    double max=-MSQ_DBL_MAX;
    size_t count = num_elements();
    for (size_t i=0; i<count; ++i) {
      Matrix3D A[MSQ_MAX_NUM_VERT_PER_ENT];
      size_t nve = elementArray[i].corner_count();
      elementArray[i].compute_corner_matrices(*this, A, nve, err);
      MSQ_ERRZERO(err);
      for (size_t j=0; j<nve; ++j) {
        min = det(A[j]) < min ? det(A[j]) : min;
        max = det(A[j]) > max ? det(A[j]) : max;
      }
    }

    if (max <= 0) {
      MSQ_SETERR(err)("Sigma_max is not positive.", MsqError::INVALID_MESH);
      return 0;
    }
      //We set delta to zero if everything in the initial mesh is valid.
      //  This causes metrics with a barrier between valid and inverted
      //  meshes to retain that barrier.  If there is a negative jacobian
      //  corner in the mesh, we set delta to a small fraction of the
      //  maximum jacobian in the mesh.
    result = (min<=MSQ_MIN) ? 0.001 * max : 0;
    computedInfos[MINMAX_SIGNED_DET3D] = result;
    note_have_info(MINMAX_SIGNED_DET3D);
  }
  
  return result;
}

/*
double PatchData::get_average_Lambda_3d( MsqError &err)
{
  double avg;
  if (have_computed_info(AVERAGE_DET3D))
  {
    avg = computedInfos[AVERAGE_DET3D];
  }
  else 
  {
    avg =0.;
    int total_num_corners =0;
    Matrix3D A[MSQ_MAX_NUM_VERT_PER_ENT];
    for (size_t i=0; i<elementArray.size(); ++i) {
      int nve = elementArray[i].corner_count();
      elementArray[i].compute_corner_matrices(*this, A, nve, err); 
      MSQ_ERRZERO(err);
      total_num_corners += nve;
      for (int c=0; c<nve; ++c) {
        avg += TargetCalculator::compute_Lambda(A[c], err); 
        MSQ_ERRZERO(err);
      }
    }

    avg = avg / total_num_corners;
    computedInfos[AVERAGE_DET3D] = avg;
    note_have_info(AVERAGE_DET3D);
  }
  return avg;
}
*/


/*! \fn PatchData::reorder()
   Physically reorder the vertices and elements in the PatchData to improve
   the locality of reference.  This method implements a Reverse Breadth First 
   Search order starting with the vertex furthest from the origin.  Other
   orderings can also be implemented.
*/
void PatchData::reorder()
{
  size_t i, j;
  const size_t num_vertex = num_nodes();
  const size_t num_elem = num_elements();

  // Step 1: Clear any cached data that will be invalidated by this
  vertexNormalPointers.clear();
  normalData.clear();
  vertexDomainDOF.clear();

  // Step 2: Make sure we have vertex-to-element adjacencies
  if (!vertAdjacencyArray.size())
    generate_vertex_to_element_data();
  
  // Step 3: Do breadth-first search
  std::vector<bool> visited( num_vertex, false );
  std::vector<size_t> vertex_order( num_vertex );
  std::vector<size_t>::iterator q1_beg, q1_end, q2_end;
  q1_beg = q1_end = q2_end = vertex_order.begin();
  // Outer loop will be done once for each disconnected chunk of mesh.
  while (q1_beg != vertex_order.end())
  {
      // Find vertex furthest from the origin
    double max = -1.0;
    size_t vtx_idx = num_vertex;
    for (i = 0; i < num_vertex; ++i)
      if (!visited[i]) 
      {
        double dist = vertexArray[i].length_squared();
        if (dist > max)
        {
          max = dist;
          vtx_idx = i;
        }
      }
    assert( vtx_idx < num_vertex);
    
    *q2_end++ = vtx_idx;;
    visited[vtx_idx] = true;
    do {
      q1_end = q2_end;
      for ( ; q1_beg != q1_end; ++q1_beg)
      {
        size_t vtx_adj_offset = vertAdjacencyOffsets[*q1_beg];
        size_t vtx_adj_end    = vertAdjacencyOffsets[*q1_beg + 1];
        for (i = vtx_adj_offset; i < vtx_adj_end; ++i)
        {
          size_t elem = vertAdjacencyArray[i];
          size_t num_elem_verts = elementArray[elem].node_count();
          size_t* elem_verts = elementArray[elem].get_vertex_index_array();
          for (j = 0; j < num_elem_verts; ++j)
          {
            size_t elem_vert = elem_verts[j];
            if (!visited[elem_vert])
            {
              *q2_end++ = elem_vert;;
              visited[elem_vert] = true;
            }
          }
        }
      }
    } while (q2_end != q1_end);
  }
  
    // Step 4: vertex_order contains the list of current vertex indices
    //         in the opposite of the order that they will occur in the 
    //         reorderd patch.  The following code will construct veretx_map
    //         from vertex_order with the following properties
    //         - vertex_map will be indexed by the current vertex index and
    //           contain the new index of that vertex (inverse of vertex_order)
    //         - the vertices will be grouped by their free/slave/fixed flag.
  std::vector<size_t> vertex_map( num_vertex );
  const size_t fixed_vtx_offset = numFreeVertices + numSlaveVertices;
  size_t free_idx = 0, slave_idx = numFreeVertices, fixed_idx = fixed_vtx_offset;
  for (i = 1; i <= num_vertex; ++i)
  {
    size_t vtx_idx = vertex_order[num_vertex - i];
    if (vtx_idx < numFreeVertices) 
      vertex_map[vtx_idx] = free_idx++;
    else if(vtx_idx < fixed_vtx_offset)
      vertex_map[vtx_idx] = slave_idx++;
    else
      vertex_map[vtx_idx] = fixed_idx++;
  }
    // make sure everything adds up
  assert( free_idx == numFreeVertices );
  assert( slave_idx == fixed_vtx_offset );
  assert( fixed_idx == num_vertex );


    // Step 5: compute element permutation
    // initialize all to "num_elem" to indicate unvisited
  std::vector<size_t> element_map( num_elem, num_elem );  
  size_t elem_idx = 0;
  for (i = 1; i <= num_vertex; ++i) 
  {
    size_t vtx_idx = vertex_order[num_vertex - i];
    size_t vtx_adj_offset = vertAdjacencyOffsets[vtx_idx];
    size_t vtx_adj_end    = vertAdjacencyOffsets[vtx_idx + 1];
    for (j = vtx_adj_offset; j < vtx_adj_end; ++j)
    {
      size_t elem = vertAdjacencyArray[j];
      if (element_map[elem] == num_elem)
        element_map[elem] = elem_idx++;
    }
  }
    // make sure everything adds up
  assert( elem_idx == num_elem );

  // Step 6:  Permute the vertices
  std::vector<MsqVertex> new_vertex_array(num_vertex);
  std::vector<Mesh::VertexHandle> new_vtx_handle_array(num_vertex);
  for (i = 0; i < num_vertex; ++i) {
    size_t new_idx = vertex_map[i];
    new_vertex_array[new_idx] = vertexArray[i];
    new_vtx_handle_array[new_idx] = vertexHandlesArray[i];
  }
  vertexArray.swap(new_vertex_array);
  vertexHandlesArray.swap(new_vtx_handle_array);

  // Step 7: Permute the elements and vertex indices for the elements
  std::vector<MsqMeshEntity> new_elem_array(num_elem);
  std::vector<Mesh::ElementHandle> new_elem_handle_array(num_elem);
  for (i = 0; i < num_elem; ++i) {
    size_t  vert_count = elementArray[i].vertex_count();
    size_t* conn_array = elementArray[i].get_vertex_index_array();
    for (j = 0; j < vert_count; ++j) {
      conn_array[j] = vertex_map[conn_array[j]];
    }

    size_t new_idx = element_map[i];
    new_elem_array[new_idx] = elementArray[i];
    new_elem_handle_array[new_idx] = elementHandlesArray[i];
  }
  elementArray.swap( new_elem_array );
  elementHandlesArray.swap( new_elem_handle_array );

  // Step 8: Clear no-longer-valid vertex-to-element adjacency info.
  if (vertAdjacencyOffsets.size()) {
    vertAdjacencyOffsets.clear();
    vertAdjacencyArray.clear();
    generate_vertex_to_element_data();
  }
  
  notify_new_patch( ); 
}


/*! 
   PatchData::move_free_vertices_constrained() moves the free vertices
   (see MsqVertex::is_free() ) as specified by the search direction (dk)
   and scale factor (step_size). After being moved, the vertices are
   projected onto the appropriate geometry.  Fixed vertices are not moved
   regardless of their corresponding dk direction.
   It is often useful to use the create_coords_momento() function before
   calling this function.
   Compile with -DMSQ_DBG3 to check that fixed vertices
   are not moved with that call.

   \param dk: must be a [nb_vtx] array of Vector3D that contains
   the direction in which to move each vertex. Fixed vertices moving
   direction should be zero, although fixed vertices will not be
   moved regardless of their corresponding dk value.
   \param nb_vtx is the number of vertices to move. must corresponds
   to the number of vertices in the PatchData.
   \param step_size will multiply the moving direction given in dk
   for each vertex.
  */
void PatchData::move_free_vertices_constrained(Vector3D dk[], size_t nb_vtx,
                                               double step_size, MsqError &err)
{
  if (nb_vtx != num_free_vertices())
  {
    MSQ_SETERR(err)("The directional vector must be of length numVertices.",
                    MsqError::INVALID_ARG);
    return;
  }
  
  size_t i;
  for (i = 0; i < num_free_vertices(); ++i) 
  {
    vertexArray[i] += (step_size * dk[i]);
    snap_vertex_to_domain(i, err);
    MSQ_ERRRTN(err);
  }
}


/*! set_free_vertices_constrained is similar to 
PatchData::move_free_vertices_constrained() except the original vertex positions
are those stored in the PatchDataVerticesMemento instead of the actual vertex
position stored in the PatchData Vertex array.  The final location is stored
in the PatchData; the PatchDataVerticesMemento is unchanged.

   \param dk: must be a [nb_vtx] array of Vector3D that contains
   the direction in which to move each vertex. Fixed vertices moving
   direction should be zero, although fixed vertices will not be
   moved regardless of their corresponding dk value.
   \param nb_vtx is the number of vertices to move. must corresponds
   to the number of vertices in the PatchData.
   \param step_size will multiply the moving direction given in dk
   for each vertex.
  */
void PatchData::set_free_vertices_constrained(PatchDataVerticesMemento* memento,
                                              Vector3D dk[],
                                              size_t nb_vtx,
                                              double step_size,
                                              MsqError &err)
{
  if (nb_vtx != memento->numVertices)
  {
    MSQ_SETERR(err)(MsqError::INVALID_ARG);
    return;
  }
  
  size_t i;
  for (i = 0; i < num_free_vertices(); ++i) 
  {
    vertexArray[i] = memento->vertices[i] + (step_size * dk[i]);
    snap_vertex_to_domain(i, err);
    MSQ_ERRRTN(err);
  }
  
    // Checks that moving direction is zero for fixed vertices.
  if (MSQ_DBG(3)) {
  for (i = 0; i < num_nodes(); ++i) {
    if (dk[i].length_squared() != 0.0)
    {
      MSQ_DBGOUT(3) << "dk["<<i<<"]: " << dk[i] << endl;
      MSQ_DBGOUT(3) << "moving a fixed vertex." << endl;
    }
  }     
  }
}


/*! Finds the maximum movement (in the distance norm) of the vertices in a
  patch.  The previous vertex positions are givena as a
  PatchDataVerticesMemento (memento).  The distance squared which each
  vertex has moved is then calculated, and the largest of those distances
  is returned.  This function only considers the movement of vertices
  that are currently 'free'.
  \param memento  a memento of this patch's vertex position at some
  (prior) time in the optimization.  
      */
double PatchData::get_max_vertex_movement_squared(PatchDataVerticesMemento*
                                                  memento,
                                                  MsqError &err)
{
  double max_dist = 0.0;
  for (size_t i = 0; i < num_free_vertices(); ++i) {
    double temp_dist = (vertexArray[i] - memento->vertices[i]).length_squared();
    if (temp_dist > max_dist)
      max_dist = temp_dist;
  }
  return max_dist;
}

/*!
 */
void PatchData::set_all_vertices_soft_fixed(MsqError &/*err*/)
{
  for(size_t i=0;i<num_free_vertices();++i)
    vertexArray[i].set_soft_fixed_flag();
}

/*!
 */
void PatchData::set_free_vertices_soft_fixed(MsqError &/*err*/)
{
  for(size_t i=0;i<num_free_vertices();++i){
    if(vertexArray[i].is_free_vertex())
      vertexArray[i].set_soft_fixed_flag();
  }
}

/*!
 */
void PatchData::set_all_vertices_soft_free(MsqError &/*err*/)
  {
    for(size_t i=0;i<num_free_vertices();++i)
      vertexArray[i].remove_soft_fixed_flag();
  }
  
/*! Get coordinates of element vertices, in canonical order.

    \param elem_index The element index in the Patch
    \param coords This vector will have the coordinates appended to it.
    If necessary, make sure to clear the vector before calling the function.
  */
void PatchData::get_element_vertex_coordinates(
  size_t elem_index,
  msq_std::vector<Vector3D> &coords,
  MsqError& /*err*/)
{
    // Check index
  if (elem_index >= num_elements())
    return;
  
    // Ask the element for its vertex indices
  const size_t *vertex_indices = elementArray[elem_index].get_vertex_index_array();
  
    // Get the coords for each indicated vertex
  size_t num_verts = elementArray[elem_index].vertex_count();
  coords.reserve(coords.size() + num_verts);
  for (size_t i = 0; i < num_verts; i++)
    coords.push_back(Vector3D(vertexArray[vertex_indices[i]]));
}

/*! This is an inefficient way of retrieving vertex_indices.
    Use PatchData::get_element_array followed by 
    MsqMeshEntity::get_vertex_index_array() if you don't need
    to fill an STL vector.
*/ 
void PatchData::get_element_vertex_indices(
  size_t elem_index,
  msq_std::vector<size_t> &vertex_indices,
  MsqError& /*err*/)
{
    // Ask the element for its vertex indices
  elementArray[elem_index].get_vertex_indices(vertex_indices);
}


void PatchData::get_vertex_element_indices(size_t vertex_index,
                                           msq_std::vector<size_t> &elem_indices,
                                           MsqError &err) 
{
  size_t count, *ptr;
  ptr = get_vertex_element_adjacencies( vertex_index, count, err );
  elem_indices.resize( count );
  memcpy( &elem_indices[0], ptr, count * sizeof(size_t) );
}

void PatchData::get_vertex_element_indices(size_t vertex_index,
                                           unsigned element_dimension,
                                           msq_std::vector<size_t> &elem_indices,
                                           MsqError &err) 
{
  elem_indices.clear();
  size_t count, *ptr;
  ptr = get_vertex_element_adjacencies( vertex_index, count, err );
  for (size_t* const end = ptr + count; ptr != end; ++ptr)
  {
    const EntityTopology type = elementArray[*ptr].get_element_type();
    const unsigned dim = TopologyInfo::dimension( type );
    if (dim == element_dimension) 
      elem_indices.push_back( *ptr );
  }
}

size_t* PatchData::get_vertex_element_adjacencies( size_t vertex_index,
                                                   size_t& array_len_out,
                                                   MsqError& err )
{
    // Make sure we've got the data
  if (vertAdjacencyArray.empty())
  {
    generate_vertex_to_element_data();
  }
  
  const size_t begin = vertAdjacencyOffsets[vertex_index];
  const size_t end = vertAdjacencyOffsets[vertex_index+1];
  array_len_out = end - begin;
  return &vertAdjacencyArray[begin];
}


/*!
    \brief This function fills a vector<size_t> with the indices
    to vertices connected to the given vertex by an edge.  If vert_indices
    is not initially empty, the function will not delete the current
    contents.  Instead, it will append the new indices at the end of
    the vector.

*/
void PatchData::get_adjacent_vertex_indices(size_t vertex_index,
                                            msq_std::vector<size_t> &vert_indices,
                                            MsqError &err)
{
    //First get elems attached to vertex[vertex_index]
  vector<size_t> elem_indices;
  vector<size_t> temp_vert_indices;
  msq_std::vector<size_t>::iterator iter;
  size_t cur_vert;
  int found=0;
  get_vertex_element_indices(vertex_index, elem_indices,err);
  MSQ_ERRRTN(err);
  MsqMeshEntity* elems=get_element_array(err);
  MSQ_ERRRTN(err);
    //get nodes attached to vertex_index... with some duplication
  while(!elem_indices.empty()){
    elems[elem_indices.back()].get_connected_vertices(vertex_index, temp_vert_indices,err); 
    MSQ_ERRRTN(err);
    elem_indices.pop_back();
  }
    //eliminate duplication.
  while(!temp_vert_indices.empty()){
    cur_vert=temp_vert_indices.back();
    temp_vert_indices.pop_back();
    iter=vert_indices.begin();
    found=0;
    while(iter!=vert_indices.end() && !found){
      if(*(iter)==cur_vert)
        found=1;
      ++iter;
    }
    if(!found)
      vert_indices.push_back(cur_vert);
  }
}

/*! Fills a vector of indices into the entities array. The entities
    in the vector are connected the given entity (ent_ind) via an
    n-diminsional entity (where 'n' is a given integer).
    Thus, if n = 0, the entities must be connected via a vertex.
    If n = 1, the entities must be connected via an edge.
    If n = 2, the entities must be connected via a two-dimensional element.
    NOTE:  if n is 2 and the elements in the entity array are
    two-dimensional, no entities should meet this criterion.
    The adj_ents vector is cleared at the beginning of the call.

*/
void PatchData::get_adjacent_entities_via_n_dim(int n, size_t ent_ind,
                                                msq_std::vector<size_t> &adj_ents,
                                                MsqError &err)
{
  //reset the vector
  adj_ents.clear();
    //vertices of this entity (given by ent_ind)
  vector<size_t> verts;
    //vector to store elements attached to the vertices in verts
  vector<size_t> elem_on_vert[MSQ_MAX_NUM_VERT_PER_ENT];
    //length of above vectos
  int length_elem_on_vert[MSQ_MAX_NUM_VERT_PER_ENT];
    //get verts on this element
  get_element_vertex_indices(ent_ind, verts, err);
  int num_vert=verts.size();
  int i=0;
  int j=0;
  for(i=0;i<num_vert;++i){
      //get elements on the vertices in verts and the number of vertices
    get_vertex_element_indices(verts[i],elem_on_vert[i],err);
    MSQ_ERRRTN(err);
    length_elem_on_vert[i]=elem_on_vert[i].size();
  }
    //this_ent is the index for an entity which is a candidate to be placed
    //into adj_ents
  size_t this_ent;
    //num of times this_ent has been found in the vectors of entity indices
  int counter=0;
  i = 0;
    //loop of each vert on ent_ind
  while(i<num_vert){
      //loop over each ent connected to vert
    j=0;
    while(j<length_elem_on_vert[i]){
        //get candidate element
      this_ent=elem_on_vert[i][j];
        //if we haven't already consider this ent
      if(this_ent!=ent_ind){
          //if this_ent occurred earlier we would have already considered it
          //so start at i and j+1
        int k1=i;
        int k2=j+1;
          //this_ent has occured once so far
        counter=1;
          //while k1 < num_vert
        while(k1<num_vert){
            //loop over entries in the elem on vert vector
          while(k2<length_elem_on_vert[k1]){
              //if it matches this_ent
            if(elem_on_vert[k1][k2]==this_ent){
                //mark it as 'seen', by making it the same as ent_ind
                //i.e., the entity  passed to us.
              elem_on_vert[k1][k2]=ent_ind;
              ++counter;
                //do not look at remaining elems in this vector
              k2+=length_elem_on_vert[k1];
            }
            else
              ++k2;
          }
          ++k1;
          k2=0;
          
        }
          //if this_ent occured enough times and isn't ent_ind
        if(counter>n && this_ent!=ent_ind){
          adj_ents.push_back(this_ent);
        }
      }
      ++j;
    }
    ++i;
  }
}

    
  


/*! \fn PatchData::update_mesh(MsqError &err)

    \brief This function copies to the TSTT mesh  the changes made to the
    free vertices / elements of the PatchData object.

*/
void PatchData::update_mesh(MsqError &err)
{
  if (!myMesh) {
    MSQ_SETERR(err)("Cannot update mesh on temporary patch.\n", MsqError::INTERNAL_ERROR);
    return;
  }
  
  for (size_t i = 0; i < vertexArray.size(); ++i)
  {
    if (!vertexArray[i].is_flag_set( MsqVertex::MSQ_HARD_FIXED ))
    {
      myMesh->vertex_set_coordinates( vertexHandlesArray[i],
                                      vertexArray[i],
                                      err ); MSQ_ERRRTN(err);
    }
    myMesh->vertex_set_byte( vertexHandlesArray[i],
                             vertexArray[i].get_flags(),
                             err ); MSQ_ERRRTN(err);
  }
}

void PatchData::generate_vertex_to_element_data()
{
  MSQ_FUNCTION_TIMER( "PatchData::generate_vertex_to_element_data" );
  
    // Skip if data already exists
  if (!vertAdjacencyArray.empty())
    return;
  
    // Allocate offset array
  vertAdjacencyOffsets.resize( num_nodes() + 1 );
  memset( &vertAdjacencyOffsets[0], 0, sizeof(size_t)*vertAdjacencyOffsets.size() );
  
    // Temporarily use offsets array to hold per-vertex element count
  msq_std::vector<MsqMeshEntity>::iterator elem_iter;
  const msq_std::vector<MsqMeshEntity>::iterator elem_end = elementArray.end();
  for (elem_iter = elementArray.begin(); elem_iter != elem_end; ++elem_iter)
  {
    size_t* conn_iter = elem_iter->get_vertex_index_array();
    const size_t* conn_end = conn_iter + elem_iter->node_count();
    for ( ; conn_iter != conn_end; ++conn_iter )
      ++vertAdjacencyOffsets[*conn_iter];
  }
  
    // Convert counts to end indices.
    // When done, vertAdjacencyOffsets will contain, for each vertex,
    // one more than the *last* index for that vertex's data in the
    // adjacency array.  This is *not* the final state for this data.
    // See comments for next loop.
  msq_std::vector<size_t>::iterator off_iter = vertAdjacencyOffsets.begin();
  const msq_std::vector<size_t>::iterator off_end = vertAdjacencyOffsets.end();
  size_t prev = *off_iter;
  ++off_iter;
  for ( ; off_iter != off_end; ++off_iter)
  {
    prev += *off_iter;
    *off_iter = prev;
  }
  
    // Allocate space for element numbers
  const size_t num_vert_uses = vertAdjacencyOffsets[num_nodes()-1];
  assert( num_vert_uses == elemConnectivityArray.size() );
  vertAdjacencyArray.resize( num_vert_uses );
  
    // Fill vertAdjacencyArray, using the indices in vertAdjacencyOffsets
    // as the location to insert the next element number in
    // vertAdjacencyArray.  When done, vertAdjacenyOffsets will contain
    // the start index for each vertex, rather than one past the last
    // index.
  for (size_t i = 0; i < elementArray.size(); ++i)
  {
    size_t* conn_iter = elementArray[i].get_vertex_index_array();
    const size_t* conn_end = conn_iter + elementArray[i].node_count();
    for ( ; conn_iter != conn_end; ++conn_iter )
    {
      const size_t array_index = --vertAdjacencyOffsets[*conn_iter];
      vertAdjacencyArray[array_index] = i;
    }
  }
  
    // Last entry should be number of vertex uses (one past the
    // last index of the last vertex.)
  vertAdjacencyOffsets[num_nodes()] = num_vert_uses;
}

void PatchData::get_subpatch(size_t center_vertex_index,
                             unsigned num_adj_elem_layers,
                             PatchData &subpatch,
                             MsqError &err)
{
  unsigned i;

    // Make sure we're in range
  if (center_vertex_index >= num_free_vertices())
  {
    MSQ_SETERR(err)("Invalid index for center vertex",MsqError::INVALID_ARG);
    return;
  }
  
    // Notify any observers of the existing subpatch that the mesh 
    // in the patch is to be changed.
  subpatch.notify_new_patch( ); 

    // Get list of vertices and elements in subpatch.
    // Ultimately, end up with arrays of unique, sorted indices.
    // It is important that the vertex indices be sorted so later
    // a reverse lookup can be done using a binary search (std::lower_bound).
  msq_std::vector<size_t> elements, vertices, offsets;
  vertices.push_back( center_vertex_index );
  for (i = 0; i < num_adj_elem_layers; ++i)
  {
    elements.clear();
    for (unsigned v = 0; v < vertices.size(); ++v)
    {
      size_t num_elem;
      const size_t* vert_elems = get_vertex_element_adjacencies( vertices[v], num_elem, err );
      MSQ_ERRRTN(err);
      elements.insert( elements.end(), vert_elems, vert_elems + num_elem );
    }
    msq_std::sort( elements.begin(), elements.end() );
    elements.erase( msq_std::unique( elements.begin(), elements.end() ), elements.end() );
    
    vertices.clear();
    for (unsigned e = 0; e < elements.size(); ++e)
    {
      MsqMeshEntity& elem = element_by_index( elements[e] );
      size_t num_vert = elem.node_count();
      const size_t* elem_verts = elem.get_vertex_index_array();
      vertices.insert( vertices.end(), elem_verts, elem_verts + num_vert );
    }
    msq_std::sort( vertices.begin(), vertices.end() );
    vertices.erase( msq_std::unique( vertices.begin(), vertices.end() ), vertices.end() );
  }
  
    // Allocate space for element connectivity info.
  size_t num_vert_uses = 0;
  for (i = 0; i < elements.size(); ++i)
    num_vert_uses += element_by_index( elements[i] ).node_count();
  subpatch.elementArray.resize( elements.size() );
  subpatch.elementHandlesArray.resize( elements.size() );
  subpatch.elemConnectivityArray.resize( num_vert_uses );
  offsets.resize( elements.size() + 1 );
  
    // Construct element connectivity data in new patch,
    // and copy element type into new patch
  size_t curr_offset = 0;
  for (i = 0; i < elements.size(); ++i)
  {
    MsqMeshEntity& elem = element_by_index( elements[i] );
    subpatch.elementArray[i].set_element_type( elem.get_element_type() );
    subpatch.elementHandlesArray[i] = elementHandlesArray[elements[i]];
    const size_t* verts = elem.get_vertex_index_array();
    offsets[i] = curr_offset;
    for (unsigned j = 0; j < elem.node_count(); ++j)
    {
      subpatch.elemConnectivityArray[curr_offset++] = 
        msq_std::lower_bound( vertices.begin(), vertices.end(), verts[j] )
        - vertices.begin();
    }
  }
  offsets[i] = curr_offset;
  
    // Store index in this patch in vertex handle array of subpatch
    // so we can determine how vertices were reordered when setting
    // vertex coordinates.
  assert(sizeof(size_t) == sizeof(void*));
  subpatch.vertexHandlesArray.resize( vertices.size() );
  size_t* vert_handles = reinterpret_cast<size_t*>(&subpatch.vertexHandlesArray[0]);
  msq_std::copy( vertices.begin(), vertices.end(), vert_handles );
  
    // All vertices except vertex at center_vertex_index are fixed.
  subpatch.byteArray.resize( vertices.size() );
  for (size_t i = 0; i < vertices.size(); ++i) {
    if (vertices[i] == center_vertex_index)
      subpatch.byteArray[i] = vertexArray[vertices[i]].get_flags() | MsqVertex::MSQ_PATCH_VTX;
    else
      subpatch.byteArray[i] = (vertexArray[vertices[i]].get_flags() & ~MsqVertex::MSQ_PATCH_VTX) | MsqVertex::MSQ_HARD_FIXED;
  }
  
    // Re-order vertices and initialize other data in subpatch
  subpatch.initialize_data( &offsets[0], &subpatch.byteArray[0], err ); 
  MSQ_ERRRTN(err);
  
    // Copy vertex data into subpatch.  subpatch.vertexHandlesArray contains
    // the indices into this PatchData for each vertex, as reordered by the
    // call to initialize_data.
  subpatch.vertexArray.resize( vertices.size() );
  for (i = 0; i < vertices.size(); ++i)
  {
    size_t vert_index = (size_t)(subpatch.vertexHandlesArray[i]);
    vertices[i] = vert_index;
    subpatch.vertexHandlesArray[i] = vertexHandlesArray[vert_index];
    subpatch.vertexArray[i] = vertexArray[vert_index];
    subpatch.vertexArray[i].flags() = subpatch.byteArray[i];
  }

  subpatch.myMesh = myMesh;
  subpatch.myDomain = myDomain;
  subpatch.mappingFuncs = mappingFuncs;
  
  notify_sub_patch( subpatch, &vertices[0], &elements[0], err ); MSQ_CHKERR(err);
}

//! Adjust the position of the specified vertex so that it
//! lies on its constraining domain.  The actual domain constraint
//! is managed by the TSTT mesh implementation
void PatchData::snap_vertex_to_domain(size_t vertex_index, MsqError &err)
{
  if (domain_set())
  {
      // if not doing normal caching
    if (normalData.empty())
    {
      get_domain()->snap_to( vertexHandlesArray[vertex_index],
                             vertexArray[vertex_index] );
    }
      // entire domain is 2-D (all vertices have a single normal)
    else if (vertexNormalPointers.empty())
    {
      get_domain()->closest_point( vertexHandlesArray[vertex_index],
                                   Vector3D(vertexArray[vertex_index]),
                                   vertexArray[vertex_index],
                                   normalData[vertex_index],
                                   err ); MSQ_ERRRTN(err);
    }
      // vertex has a singel normal
    else if (2 == vertexDomainDOF[vertex_index])
    {
      get_domain()->closest_point( vertexHandlesArray[vertex_index],
                                   Vector3D(vertexArray[vertex_index]),
                                   vertexArray[vertex_index],
                                   *vertexNormalPointers[vertex_index],
                                   err ); MSQ_ERRRTN(err);
    }
      // vertex as multiple normals cached, one per 2-D element
    else
    {
        // get closest point
      get_domain()->snap_to( vertexHandlesArray[vertex_index],
                             vertexArray[vertex_index] );

        // update per-element normal values
      offsetArray.clear(); // use this array for tmp storage 
      get_vertex_element_indices( vertex_index, 2, offsetArray, err );
      MSQ_ERRRTN(err);
      if (offsetArray.size())
      {
        for (size_t i = 0; i < offsetArray.size(); ++i)
        {
          vertexNormalPointers[vertex_index][i] = vertexArray[vertex_index];
          offsetArray[i] = (size_t)(elementHandlesArray[offsetArray[i]]);
        }
        assert( sizeof(size_t) == sizeof(Mesh::ElementHandle) );
        get_domain()->normal_at( (Mesh::ElementHandle*)(&offsetArray[0]), 
                                 vertexNormalPointers[vertex_index],
                                 offsetArray.size(), err );
        MSQ_ERRRTN(err);
      }
    }
  }
}


void PatchData::update_cached_normals( MsqError& err )
{
  size_t i;
  
  MeshDomain* domain = get_domain();
  if (!domain)
  {
    MSQ_SETERR(err)( "No domain constraint set.", MsqError::INVALID_STATE );
    return;
  }
  
    // Determine which vertices lie on surfaces
  vertexDomainDOF.resize( num_nodes() );
  domain->domain_DoF( &vertexHandlesArray[0], &vertexDomainDOF[0], num_nodes(), err );
  MSQ_ERRRTN(err);
  
    // Count how many vertices have a single normal
  // Sun doesn't support partial template specialization, so can't use std::count
  //size_t n = msq_std::count( vertexDomainDOF.begin(), vertexDomainDOF.end(), 2 );
  size_t n = 0;
  msq_std::vector<unsigned short>::iterator k;
  for ( k = vertexDomainDOF.begin(); k != vertexDomainDOF.end(); ++k)
    if (*k == 2)
      ++n;
  
    // If all vertices are on a surface, pass in the existing handles array
    // and store a single normal per vertex.
  if (n == num_nodes())
  {
    normalData.resize( n );
    msq_std::copy( vertexArray.begin(), vertexArray.end(), normalData.begin() );
    domain->normal_at( &vertexHandlesArray[0], &normalData[0], num_nodes(), err );
    vertexNormalPointers.clear();
    MSQ_ERRRTN(err);
//NOTE: stop here if all vertices are on a surface
    return;
  }
    
    // count how many per-element vertex normals need to be stored
  msq_std::vector<size_t> vertex_elements;
  const size_t nn = num_nodes();
  size_t num_normals = n;
  for (i = 0; i < nn; ++i)
  {
    if (vertexDomainDOF[i] != 2)
    {
      get_vertex_element_indices( i, 2, vertex_elements, err ); MSQ_ERRRTN(err);
      num_normals += vertex_elements.size();
    }
  }
  normalData.resize( num_normals );
  vertexNormalPointers.reserve( num_nodes() );
  msq_std::vector<Vector3D>::iterator norm_iter = normalData.begin();
  
    // Fill a temporary handle array with entity handles and
    // the normalData array with the corresponding vertex coordinates
    // for a call to get the normals.  Also, populate vertexNormalPointers
    // array with offsets into normalData.
  msq_std::vector<Mesh::EntityHandle> handles;
  handles.reserve( num_normals );
  for (i = 0; i < nn; ++i)
  {
    vertexNormalPointers.push_back( &*norm_iter );
    if (vertexDomainDOF[i] == 2)
    {
      *norm_iter = vertexArray[i];
      handles.push_back( vertexHandlesArray[i] );
      ++norm_iter;
    }
    else
    {
      get_vertex_element_indices( i, 2, vertex_elements, err ); MSQ_ERRRTN(err);
      for (msq_std::vector<size_t>::iterator j = vertex_elements.begin();
           j != vertex_elements.end(); ++j)
      {
        handles.push_back( elementHandlesArray[*j] );
        *norm_iter = vertexArray[i];
        ++norm_iter;
      }
    }
  }
  assert( norm_iter == normalData.end() );
  
    // Get the actual normals
  domain->normal_at( &handles[0], &normalData[0], num_normals, err ); MSQ_ERRRTN(err);
}


void PatchData::get_domain_normal_at_element(size_t elem_index,
                                             Vector3D &surf_norm,
                                             MsqError &err) const
{
  if (domain_set())
  {
    elementArray[elem_index].get_centroid(surf_norm, *this, err); MSQ_ERRRTN(err);
    get_domain()->normal_at( elementHandlesArray[elem_index], surf_norm );
  }
  else
    MSQ_SETERR(err)( "No domain constraint set.", MsqError::INVALID_STATE );
}


void PatchData::get_domain_normal_at_mid_edge( size_t elem_index,
                                               unsigned edge_num,
                                               Vector3D& normal,
                                               MsqError& err ) const
{
  if (!domain_set()) {
    MSQ_SETERR(err)("No domain constraint set.", MsqError::INVALID_STATE );
    return;
  }
  
  const MsqMeshEntity& elem = element_by_index( elem_index );
  const unsigned* edge = TopologyInfo::edge_vertices( elem.get_element_type(),
                                                      edge_num, err );
                                                      MSQ_ERRRTN(err);
  const MsqVertex& v1 = vertex_by_index( elem.get_vertex_index_array()[edge[0]] );
  const MsqVertex& v2 = vertex_by_index( elem.get_vertex_index_array()[edge[1]] );
  normal = 0.5 * (v1 + v2);
  get_domain()->normal_at( elementHandlesArray[elem_index], normal );
}

void PatchData::get_domain_normals_at_corners( size_t elem_index,
                                               Vector3D normals_out[],
                                               MsqError& err ) 
{
  if (!domain_set())
  {
    MSQ_SETERR(err)( "No domain constraint set.", MsqError::INVALID_STATE );
    return;
  }
  
  if (2 != TopologyInfo::dimension( elementArray[elem_index].get_element_type() ))
  {
    MSQ_SETERR(err)( "Attempt to get corners of non-surface element", MsqError::INVALID_ARG );
    return;
  }
  
  if (normalData.empty())
  {
    update_cached_normals( err ); MSQ_ERRRTN(err);
  }
  
  MsqMeshEntity& elem = elementArray[elem_index];
  const unsigned count = elem.vertex_count();
  const size_t* const vertex_indices = elem.get_vertex_index_array();
  for (size_t i = 0; i < count; ++i)
  {
    const size_t v = vertex_indices[i];
    if (vertexNormalPointers.empty())
    {
      normals_out[i] = normalData[v];
    }
    else if(vertexDomainDOF[v] == 2)
    {
      normals_out[i] = *vertexNormalPointers[v];
    }
    else
    {
      offsetArray.clear();  // use this temp array instead of allocating a new one.
      get_vertex_element_indices( v, 2, offsetArray, err ); MSQ_ERRRTN(err);
      
      msq_std::vector<size_t>::iterator iter;
      iter = msq_std::find( offsetArray.begin(), offsetArray.end(), elem_index );
      if (iter == offsetArray.end())
      {
        MSQ_SETERR(err)("Invalid element handle (not 2D?)", MsqError::INVALID_ARG );
        return;
      }
      normals_out[i] = vertexNormalPointers[v][iter - offsetArray.begin()];
      offsetArray.clear();
    }
  }
}  

void PatchData::get_domain_normal_at_corner( size_t elem_index,
                                             unsigned corner,
                                             Vector3D& normal,
                                             MsqError& err ) 
{
  if (!domain_set())
  {
    MSQ_SETERR(err)( "No domain constraint set.", MsqError::INVALID_STATE );
    return;
  }
  
  if (2 != TopologyInfo::dimension( elementArray[elem_index].get_element_type() ))
  {
    MSQ_SETERR(err)( "Attempt to get corners of non-surface element", MsqError::INVALID_ARG );
    return;
  }
  
  if (normalData.empty())
  {
    update_cached_normals( err ); MSQ_ERRRTN(err);
  }
  
  MsqMeshEntity& elem = elementArray[elem_index];
  const size_t* const vertex_indices = elem.get_vertex_index_array();
  const size_t v = vertex_indices[corner];
  if (vertexNormalPointers.empty())
  {
    normal = normalData[v];
  }
  else if(vertexDomainDOF[v] == 2)
  {
    normal = *vertexNormalPointers[v];
  }
  else
  {
    offsetArray.clear();  // use this temp array instead of allocating a new one.
    get_vertex_element_indices( v, 2, offsetArray, err ); MSQ_ERRRTN(err);

    msq_std::vector<size_t>::iterator iter;
    iter = msq_std::find( offsetArray.begin(), offsetArray.end(), elem_index );
    if (iter == offsetArray.end())
    {
      MSQ_SETERR(err)("Invalid element handle (not 2D?)", MsqError::INVALID_ARG );
      return;
    }
    normal = vertexNormalPointers[v][iter - offsetArray.begin()];
    offsetArray.clear();
  }
}  
  

void PatchData::set_mesh(Mesh* ms)        
{ 
  myMesh = ms; 
    // observers should treat this the same as if the
    // instance of this object wzs being deleted.
  notify_patch_destroyed();
}

void PatchData::set_domain(MeshDomain* d) 
{ 
  myDomain = d;
  
    // clear all cached data from the previous domain
  vertexNormalPointers.clear();
  normalData.clear();
  vertexDomainDOF.clear();

    // observers should treat this the same as if the
    // instance of this object wzs being deleted.
  notify_patch_destroyed();
}

ostream& operator<<( ostream& stream, const PatchData& pd )
{
   size_t i;
   
   stream << "Vertices: " << endl;
   for (i = 0; i < pd.num_nodes(); ++i)
   {
      stream << i << ". " 
             << pd.vertexHandlesArray[i] << " (" 
             << pd.vertexArray[i].x() << ","
             << pd.vertexArray[i].y() << ","
             << pd.vertexArray[i].z()
             << ") ";
      if (pd.vertexArray[i].is_flag_set( MsqVertex::MSQ_CULLED ))
        stream << "S";
      if (pd.vertexArray[i].is_flag_set( MsqVertex::MSQ_HARD_FIXED ))
        stream << "H";
      
      if (pd.vertAdjacencyArray.size())
      {
        size_t j = pd.vertAdjacencyOffsets[i];
        size_t end = pd.vertAdjacencyOffsets[i+1];
        for ( ; j < end; ++j )
          stream << " " << pd.vertAdjacencyArray[j];
      }
      
      stream << endl;
   }
   
   stream << "Elements: " << endl;
   for (i = 0; i < pd.num_elements(); ++i)
   {
      stream << i << ". " << pd.elementHandlesArray[i] << " ";
      switch (pd.elementArray[i].get_element_type()) {
        case POLYGON:       stream << "Polygon";    break;
        case TRIANGLE:      stream << "Tri";        break;
        case QUADRILATERAL: stream << "Quad";       break;
        case POLYHEDRON:    stream << "Polyhedron"; break;
        case TETRAHEDRON:   stream << "Tet";        break;
        case HEXAHEDRON:    stream << "Hex";        break;
        case PRISM:         stream << "Wedge";      break;
        case PYRAMID:       stream << "Pyr";        break;
        default:            stream << "Unknown";    break;
      }
      stream << pd.elementArray[i].node_count() << ": ";
      for (size_t j = 0; j < pd.elementArray[i].node_count(); ++j)
        stream << pd.elementArray[i].get_vertex_index_array()[j] << " ";
      stream << endl;
   }
   stream << endl;
   
   stream << "Mesh: " << (pd.myMesh?"yes":"no") << endl;
   stream << "Domain: " << (pd.myDomain?"yes":"no") << endl;
//   stream << "mType: " << (pd.mType==PatchData::VERTICES_ON_VERTEX_PATCH?"vert-on-vert":
//                           pd.mType==PatchData::ELEMENTS_ON_VERTEX_PATCH?"elem-on-vert":
//                           pd.mType==PatchData::GLOBAL_PATCH?"global":"unknown") << endl;
   
   if (pd.haveComputedInfos)
   {
     stream << "ComputedInfos:" << endl;
     if (pd.have_computed_info(PatchData::MIN_UNSIGNED_AREA))
       stream << "\t MIN_UNSINGED_AREA = " << pd.computedInfos[PatchData::MIN_UNSIGNED_AREA] << endl;
     if (pd.have_computed_info(PatchData::MAX_UNSIGNED_AREA))
       stream << "\t MAX_UNSIGNED_AREA = " << pd.computedInfos[PatchData::MAX_UNSIGNED_AREA] << endl;
     if (pd.have_computed_info(PatchData::MIN_EDGE_LENGTH))
       stream << "\t MIN_EDGE_LENGTH = " << pd.computedInfos[PatchData::MIN_EDGE_LENGTH] << endl;
     if (pd.have_computed_info(PatchData::MAX_EDGE_LENGTH))
       stream << "\t MAX_EDGE_LENGTH = " << pd.computedInfos[PatchData::MAX_EDGE_LENGTH] << endl;
     if (pd.have_computed_info(PatchData::MINMAX_SIGNED_DET2D))
       stream << "\t MINMAX_SIGNED_DET2D = " << pd.computedInfos[PatchData::MINMAX_SIGNED_DET2D] << endl;
     if (pd.have_computed_info(PatchData::MINMAX_SIGNED_DET3D))
       stream << "\t MINMAX_SIGNED_DET3D = " << pd.computedInfos[PatchData::MINMAX_SIGNED_DET3D] << endl;
     if (pd.have_computed_info(PatchData::AVERAGE_DET3D))
       stream << "\t AVERAGE_DET3D = " << pd.computedInfos[PatchData::AVERAGE_DET3D] << endl;
  }
  
  return stream << endl;
}

void print_patch_data( const PatchData& pd )
{
  std::cout << pd << std::endl;
}

void PatchData::initialize_patch( EntityTopology* elem_type_array,
                                  size_t* elem_offset_array,
                                  unsigned char* vertex_flags,
                                  MsqError& err )
{
  size_t i;
  
    // set element types
  elementArray.resize( elementHandlesArray.size() );
  for (i = 0; i < elementHandlesArray.size(); ++i)
    elementArray[i].set_element_type( elem_type_array[i] );
  
    // get element connectivity, group vertices by free/slave/fixed state
  initialize_data( elem_offset_array, vertex_flags, err ); MSQ_ERRRTN(err);
  
    // get vertex coordinates
  vertexArray.resize( vertexHandlesArray.size() );
  get_mesh()->vertices_get_coordinates( &vertexHandlesArray[0],
                                    &vertexArray[0],
                                    vertexHandlesArray.size(),
                                    err ); MSQ_ERRRTN(err);

    // set vertex flags
  for (i = 0; i < vertexArray.size(); ++i)
    vertexArray[i].flags() = vertex_flags[i];

    // get culled bit from mesh
  get_mesh()->vertices_get_byte( &vertexHandlesArray[0],
                                 vertex_flags,
                                 vertexHandlesArray.size(),
                                 err ); MSQ_ERRRTN(err);

  for (i = 0; i < vertexArray.size(); ++i)
    vertexArray[i].flags() |= (vertex_flags[i] & MsqVertex::MSQ_CULLED);
}

void PatchData::initialize_data( size_t* elem_offset_array, 
                                 unsigned char* vertex_flags,
                                 MsqError& err )
{
    // Clear out data specific to patch
  vertexNormalPointers.clear();
  normalData.clear();
  vertexDomainDOF.clear();
  
    // Clear any vertex->element adjacency data.  It
    // is probably invalid, and certainly will be by the time
    // this function completes if the mesh contains higher-order
    // elements.
  vertAdjacencyArray.clear();
  vertAdjacencyOffsets.clear();
  size_t i, j;
  bool have_slave_vertices = false;
  for (i = 0; i < elementArray.size(); ++i)
  {
    size_t start = elem_offset_array[i];
    size_t conn_len = elem_offset_array[i+1] - start;
    assert(conn_len > 0);
    elementArray[i].set_connectivity( &elemConnectivityArray[start], conn_len );
    if (conn_len > elementArray[i].vertex_count())
      have_slave_vertices = true;
  }
    
    // Use vertAdjacencyOffsets array as temporary storage.
  vertAdjacencyOffsets.resize( vertexHandlesArray.size() + 1 );
  size_t* vertex_index_map = &vertAdjacencyOffsets[0];
  
    // Count number of free vertices and initialize vertex_index_map
  numFreeVertices = 0;
  for (i = 0; i < vertexHandlesArray.size(); ++i) {
    if (!(vertex_flags[i] & MsqVertex::MSQ_HARD_FIXED))
      ++numFreeVertices;
    vertex_index_map[i] = i;
  }
  
    // Re-order vertices such that all free vertices are
    // first in the list.  Construct map from old to new
    // position in list for updating element connectivity.
  i = 0;
  j = numFreeVertices;
  for (;; ++i, ++j) {
    for (; i < numFreeVertices && !(vertex_flags[i] & MsqVertex::MSQ_HARD_FIXED); ++i);
    if (i == numFreeVertices)
      break;
    for ( ; vertex_flags[j] & MsqVertex::MSQ_HARD_FIXED; ++j);
    msq_std::swap( vertex_index_map[i], vertex_index_map[j] );
    msq_std::swap( vertexHandlesArray[i], vertexHandlesArray[j] );
    msq_std::swap( vertex_flags[i], vertex_flags[j] );
  }
  assert( i == numFreeVertices );
  assert( j <= vertexHandlesArray.size() );
  
    // Reorder vertices such that free, slave vertices
    // occur after free, non-slave vertices in list.
    // CHANGE ME: for now, assume all higher-order nodes are slaves.
  numSlaveVertices = 0;
  if (have_slave_vertices) {
      // mark slave vertices, and count them.
    for (i = 0; i < elementArray.size(); ++i)
    {
      size_t start = elem_offset_array[i];
      size_t conn_len = elem_offset_array[i+1] - start;
      for (j = elementArray[i].vertex_count(); j < conn_len; ++j) {
        const size_t vert_idx = vertex_index_map[elemConnectivityArray[start+j]];
        assert(vert_idx < vertexHandlesArray.size());
        if (vert_idx < numFreeVertices && !(vertex_flags[i]&MsqVertex::MSQ_DEPENDENT)) {
          vertex_flags[i] |= MsqVertex::MSQ_DEPENDENT;
          ++numSlaveVertices;
        }
      }
    }
    numFreeVertices -= numSlaveVertices;

      // Re-order free vertices such that all slave vertices are
      // last in the list.  Construct map from old to new
      // position in list for updating element connectivity.
    i = 0;
    j = numFreeVertices;
    for (;; ++i, ++j) {
      for (; i < numFreeVertices && !(vertex_flags[i]&MsqVertex::MSQ_DEPENDENT); ++i);
      if (i >= numFreeVertices)
        break;
      for ( ; vertex_flags[j] & MsqVertex::MSQ_DEPENDENT; ++j);
      msq_std::swap( vertex_index_map[i], vertex_index_map[j] );
      msq_std::swap( vertexHandlesArray[i], vertexHandlesArray[j] );
    }
    assert( i == numFreeVertices );
    assert( j <= numFreeVertices + numSlaveVertices );
  }  
  
    // Update element connectivity data for new vertex indices
  for (i = 0; i < elemConnectivityArray.size(); ++i)
    elemConnectivityArray[i] = vertex_index_map[elemConnectivityArray[i]];

    // Clear temporary data    
  vertAdjacencyOffsets.clear();
  
  notify_new_patch( );
}

size_t PatchData::num_corners() const
{
  size_t result = 0;
  for (unsigned i = 0; i < elementArray.size(); ++i)
    result += elementArray[i].vertex_count();
  return result;
}


void PatchData::fill( size_t num_vertex, const double* coords,
                      size_t num_elem, EntityTopology type, 
                      const size_t* connectivity,
                      const bool* fixed, 
                      MsqError& err )
{
  msq_std::vector<EntityTopology> types(num_elem);
  msq_std::fill( types.begin(), types.end(), type );
  this->fill( num_vertex, coords, num_elem, &types[0], connectivity, fixed, err );
  MSQ_CHKERR(err);
}
  
    
void PatchData::fill( size_t num_vertex, const double* coords,
                      size_t num_elem, const EntityTopology* types,
                      const size_t* conn,
                      const bool* fixed, 
                      MsqError& err )
{
  msq_std::vector<size_t> lengths( num_elem );
  msq_std::transform( types, types + num_elem, lengths.begin(), 
                      std::ptr_fun(TopologyInfo::corners) );
  this->fill( num_vertex, coords, num_elem, types, &lengths[0], conn, fixed, err );
  MSQ_CHKERR(err);
} 
    
void PatchData::fill( size_t num_vertex, const double* coords,
                      size_t num_elem, const EntityTopology* types,
                      const size_t* lengths,
                      const size_t* conn,
                      const bool* fixed,
                      MsqError& err )
{
  size_t i;
  
    // count vertex uses
  size_t num_uses = msq_std::accumulate( lengths, lengths + num_elem, 0 );

    // Allocate storage for data
  vertexArray.resize( num_vertex );
  vertexHandlesArray.resize( num_vertex );
  elementArray.resize( num_elem );
  elementHandlesArray.resize( num_elem );
  elemConnectivityArray.resize( num_uses );
  numFreeVertices = 0;
  numSlaveVertices = 0;
  
  if (fixed) {
    byteArray.resize( num_vertex );
    for (i = 0; i < num_vertex; ++i)
      byteArray[i] = fixed[i] ? MsqVertex::MSQ_HARD_FIXED : MsqVertex::MSQ_PATCH_VTX;
  }
  else {
    byteArray.clear();
    byteArray.resize( num_vertex, MsqVertex::MSQ_PATCH_VTX );
  }
  
  for (i = 0; i < num_elem; ++i) {
    element_by_index(i).set_element_type( types[i] );
    elementHandlesArray[i] = (Mesh::ElementHandle)i;
  }
  for (i = 0; i < num_vertex; ++i)
    vertexHandlesArray[i] = (Mesh::VertexHandle)i;
  
  memcpy( get_connectivity_array(), conn, num_uses * sizeof(size_t) );
  
  msq_std::vector<size_t> offsets( num_elem + 1 );
  size_t sum = offsets[0] = 0;
  for (i = 1; i <= num_elem; ++i)
    offsets[i] = sum += lengths[i-1];
  
  this->initialize_data( &offsets[0], &byteArray[0], err );  MSQ_ERRRTN(err);
  
    // initialize_data will re-order vertex handles and
    // update element connectivity accordingly.  Use
    // the values we stored in vertexHandlesArray to
    // figure out the new index of each vertex, and initialize
    // the vertex.
  for (i = 0; i < num_vertex; ++i) 
    vertex_by_index((size_t)vertexHandlesArray[i]) = coords + 3*i;
  
  for (i = 0; i < num_vertex; ++i)
    vertex_by_index(i).flags() = byteArray[i];
}
  
void PatchData::make_handles_unique( Mesh::EntityHandle* handles,
                                     size_t& count,
                                     size_t* index_map )
{
  if (count < 2)
  {
    return;
  }
    // save this now, as we'll be changing count later
  const size_t* index_end = index_map + count;
  
  if (index_map)
  {
      // Copy input handle list into index map as a temporary
      // copy of the input list.
    assert( sizeof(Mesh::EntityHandle) == sizeof(size_t) );
    memcpy( index_map, handles, count*sizeof(size_t) );
  }

    // Make handles a unique list
  msq_std::sort( handles, handles + count );
  Mesh::EntityHandle* end = msq_std::unique( handles, handles + count );
  count = end - handles;

  if (index_map)
  {
      // Replace each handle in index_map with the index of
      // its position in the handles array
    Mesh::EntityHandle* pos;
    for (size_t* iter = index_map; iter != index_end; ++iter)
    {
      pos = msq_std::lower_bound( handles, handles + count, (Mesh::EntityHandle)*iter );
      *iter = pos - handles;
    }
  }
}

void PatchData::fill_global_patch( MsqError& err )
{
  GlobalPatch gp;
  gp.set_mesh( get_mesh() );
  PatchIterator iter( &gp );
  bool b = iter.get_next_patch( *this, err ); MSQ_ERRRTN(err);
  assert(b);
}


void PatchData::set_mesh_entities( 
                          msq_std::vector<Mesh::ElementHandle>& elements,
                          msq_std::vector<Mesh::VertexHandle>& free_vertices,
                          MsqError& err )
{
  if (!get_mesh()) {
    MSQ_SETERR(err)("No Mesh associated with PatchData.", MsqError::INVALID_STATE );
    return;
  }
  
  elementHandlesArray = elements;
  get_mesh()->elements_get_attached_vertices( &elementHandlesArray[0],
                                              elementHandlesArray.size(),
                                              vertexHandlesArray,
                                              offsetArray,
                                              err ); MSQ_ERRRTN(err);
  
    // Construct CSR-rep element connectivity data
  size_t num_vert = vertexHandlesArray.size();
  elemConnectivityArray.resize( num_vert );
  make_handles_unique( &vertexHandlesArray[0], 
                       num_vert, 
                       &elemConnectivityArray[0] ); 
  vertexHandlesArray.resize( num_vert );

    // Get element topologies
  msq_std::vector<EntityTopology> elem_topologies(elementHandlesArray.size());
  get_mesh()->elements_get_topologies( &elementHandlesArray[0],
                                       &elem_topologies[0],
                                       elementHandlesArray.size(),
                                       err ); MSQ_ERRRTN(err);
  
  byteArray.resize( vertexHandlesArray.size() );
  bool* fixed_flags;
  if (sizeof(bool) == sizeof(unsigned char)) {
      // sizeof(bool) == sizeof(unsigned char) on every platform I'm
      // aware of.  However, the standard doesn't guarantee it so this
      // code is here for portability.
    fixed_flags = reinterpret_cast<bool*>(&byteArray[0]);
  }
  else {
    fixed_flags = new bool[vertexHandlesArray.size()];
  }
  
    // Get vertex-fixed flags
  get_mesh()->vertices_get_fixed_flag( &vertexHandlesArray[0], 
                                       fixed_flags,
                                       vertexHandlesArray.size(),
                                       err);
  
  if (sizeof(bool) != sizeof(unsigned char)) {
    for (size_t i = 0; i < vertexHandlesArray.size(); ++i) 
      if (fixed_flags[i])
        byteArray[i] = MsqVertex::MSQ_HARD_FIXED;
      else
        byteArray[i] = 0;
    delete fixed_flags;
  }
  MSQ_ERRRTN(err);
  
    // if free_vertices is not empty, then we need to mark as
    // fixed any vertices *not* in that list.
  if (free_vertices.empty()) {
    for (size_t i = 0; i < vertexHandlesArray.size(); ++i) {
      assert( byteArray[i] <= 1 ); // make sure cast from bool works as expected (true == 1)
      byteArray[i] |= MsqVertex::MSQ_PATCH_VTX;
    }
  }
  else {
      // sort and remove duplicates from free_vertices list.
    msq_std::sort(free_vertices.begin(), free_vertices.end());
    free_vertices.erase( 
        msq_std::unique(free_vertices.begin(), free_vertices.end()), 
        free_vertices.end() );
    
    for (size_t i = 0; i < vertexHandlesArray.size(); ++i) {
      assert( byteArray[i] <= 1 ); // make sure cast from bool works as expected (true == 1)
      if (msq_std::binary_search(free_vertices.begin(), 
                                 free_vertices.end(), 
                                 vertexHandlesArray[i]))
        byteArray[i] |= MsqVertex::MSQ_PATCH_VTX;
      else
        byteArray[i] |= MsqVertex::MSQ_HARD_FIXED;
    }
  }

    // Store element connectivty, and do any additional
    // data arrangement that needs to be done.
  initialize_patch( &elem_topologies[0], &offsetArray[0], &byteArray[0], err ); 
  MSQ_CHKERR(err);
}


void PatchData::set_mapping_functions( MappingFunctionSet* mfs )
  { mappingFuncs = mfs; }

MappingFunctionSet* PatchData::get_mapping_functions() const
  { return mappingFuncs; }

const MappingFunction* PatchData::get_mapping_function( EntityTopology type ) const
{
  return mappingFuncs ? mappingFuncs->get_function(type) : 0;
}

void PatchData::get_sample_location( size_t element_index,
                              unsigned sample_dim,
                              unsigned sample_num,
                              Vector3D& result,
                              MsqError& err ) const
{
  const MsqMeshEntity& elem = element_by_index( element_index );
  const unsigned ho_bits = higher_order_node_bits( element_index );
  const MappingFunction* const f = get_mapping_function( elem.get_element_type() );
  if (!f) {
    MSQ_SETERR(err)("No mapping function", MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  
  std::vector<double> coeff;
  switch (sample_dim) {
    case 0:
      f->coefficients_at_corner( sample_num, ho_bits, coeff, err );
      break;
    case 1:
      f->coefficients_at_mid_edge( sample_num, ho_bits, coeff, err );
      break;
    case 2:
      if (TopologyInfo::dimension( elem.get_element_type() ) != 2) {
        f->coefficients_at_mid_face( sample_num, ho_bits, coeff, err );
        break;
      }
    case 3:
      f->coefficients_at_mid_elem( ho_bits, coeff, err );
      break;
  }
  MSQ_ERRRTN( err );
  
  const size_t* const conn = elem.get_vertex_index_array();
  result[2] = result[1] = result[0] = 0.0;
  for (unsigned i = 0; i < coeff.size(); ++i) 
    result += coeff[i] * vertex_by_index(conn[i]);
}


unsigned PatchData::higher_order_node_bits( size_t element_index ) const
{
  const MsqMeshEntity& elem = element_by_index(element_index);
  EntityTopology type = elem.get_element_type();
  const size_t* conn = elem.get_vertex_index_array();
  const size_t n = elem.node_count();

  MsqError err;
  bool have_midedge, have_midface, have_midelem;
  unsigned num_edge = 0, num_face = 0, num_corner = TopologyInfo::corners(type);
  TopologyInfo::higher_order( type, n, have_midedge, have_midface, have_midelem, err );
  num_edge = TopologyInfo::edges(type);
  if (TopologyInfo::dimension(type) == 2) {
    have_midelem = have_midface;
    have_midface = false;
  }
  else 
    num_face = TopologyInfo::faces(type);
  
  unsigned result = 0;
  if (have_midedge) {
    for (unsigned i = 0; i < num_edge; ++i) {
      if (!(vertex_by_index(conn[num_corner+i]).get_flags() & MsqVertex::MSQ_DEPENDENT))
        result |= (1 << i);
    }
  }
  if (have_midface) {
    for (unsigned i = 0; i < num_face; ++i) {
      if (!(vertex_by_index(conn[num_corner+num_edge+i]).get_flags() & MsqVertex::MSQ_DEPENDENT))
        result |= (1 << (i+num_edge));
    }
  }
  if (have_midelem && !(vertex_by_index(conn[num_corner+num_edge+num_face]).get_flags() & MsqVertex::MSQ_DEPENDENT))
    result |= (1 << (num_edge+num_face));

  return result;
}

bool PatchData::attach_extra_data( ExtraData* data )
{
  if (data->patchNext) {
    return false;
  }
  
  if (!data->patchPtr)
    data->patchPtr = this;
  else if(data->patchPtr != this)
    return false;
  
  data->patchNext = dataList;
  dataList = data;
  return true;
}

bool PatchData::remove_extra_data( ExtraData* data )
{
  if (data->patchPtr != this)
    return false;
  
  if (dataList == data) {
    dataList = data->patchNext;
    data->patchNext = 0;
    data->patchPtr = 0;
    return true;
  }
  
  for (ExtraData* iter = dataList; iter; iter = iter->patchNext)
    if (iter->patchNext == data) {
      iter->patchNext = data->patchNext;
      data->patchNext = 0;
      data->patchPtr = 0;
      return true;
    }
  
  return false;
}

void PatchData::notify_new_patch()
{
  for (ExtraData* iter = dataList; iter; iter = iter->patchNext)
    iter->notify_new_patch();
}

void PatchData::notify_sub_patch( PatchData& sub_patch,
                                  const size_t* vertex_map,
                                  const size_t* element_map,
                                  MsqError& err )
{
  for (ExtraData* iter = dataList; iter; iter = iter->patchNext)
  {
    iter->notify_sub_patch( sub_patch, vertex_map, element_map, err );
    MSQ_ERRRTN(err);
  }
}

void PatchData::notify_patch_destroyed()
{
    // Remove all ExtraDatas from list and notify them
    // that they are being destroyed.
  while( dataList ) {
    ExtraData* dead_data = dataList;
    dataList = dataList->patchNext;
    dead_data->patchNext = 0;
    dead_data->patchPtr = 0;
    dead_data->notify_patch_destroyed( );
  }
}


} // namespace Mesquite
