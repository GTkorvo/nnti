/*!
  \file   PatchData.cpp

  \author Thomas Leurent
  \author Michael Brewer
  \date   2002-01-17
*/

#include "PatchData.hpp"
#include "MsqMeshEntity.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MeshSet.hpp"

#include <list>
#include "MsqMessage.hpp"
#include "MsqTimer.hpp"
using namespace Mesquite;  
using namespace std;

#undef __FUNC__
#define __FUNC__ "PatchData::PatchData"
PatchData::PatchData()
  : meshSet(NULL),
    domainSet(false),
    mType(UNDEFINED_PATCH_TYPE),
    numVertices(0),
    numElements(0),
    vertexArray(NULL),
    vertexHandlesArray(NULL),
    elementArray(NULL),
    elementHandlesArray(NULL),
    v2E(NULL),
    v2eOffset(NULL),
    subpatchIndexArray(NULL),
    vertexArraySize(0),
    elemArraySize(0),
    v2eSize(0),
    v2eOffsetSize(0),
    v2eValid(false),
    subpatchIndexSize(0)
{}


// Destructor
#undef __FUNC__
#define __FUNC__ "PatchData::~PatchData" 
PatchData::~PatchData()
{
  delete [] vertexArray;
  delete [] vertexHandlesArray;
  delete [] elementArray;
  delete [] elementHandlesArray;
  delete [] v2E;
  delete [] v2eOffset;
  delete [] subpatchIndexArray;
}

/*! \fn PatchData::reorder()
   Physically reorder the vertices and elements in the PatchData to improve
   the locality of reference.  This method implements a Reverse Breadth First 
   Search order starting with the vertex furthest from the origin.  Other
   orderings can also be implemented.
*/
void PatchData::reorder()
{
  const size_t numv = numVertices;
  const size_t nume = numElements;

  size_t *vtx;
  size_t *tmp;

  size_t *sta = new size_t[numv + 1];
  size_t *vte;
  size_t *ord = new size_t[numv];
  size_t *per = new size_t[numv];
  size_t *pel;
  size_t *que1 = new size_t[numv];
  size_t *que2 = new size_t[numv];

  MsqVertex *v2a;
  Mesh::VertexHandle *v2h;
  MsqMeshEntity *e2a;
  Mesh::ElementHandle *e2h;

  double val, max;

  size_t toc;
  size_t vtc;
  size_t idx;
  size_t loc;
  size_t i, j;
  size_t q1l, q2l, q;
  size_t st, en;

  // Step 0:  Make sure patch data is valid.

  // Step 1:  Find the length of the element to vertex list for each 
  //          individual vertex.

  memset(sta, 0, (numv+1)*sizeof(size_t));
  for (i = 0; i < nume; ++i) {
    vtc = elementArray[i].vertex_count();
    vtx = elementArray[i].get_modifiable_vertex_index_array();

    for (j = 0; j < vtc; ++j) {
      ++sta[vtx[j]];
    }
  }

  // Step 2:  Compute the offsets, total length of the element to vertex
  //          list, and allocate the data.

  toc = sta[0];
  sta[0] = 0;
  for (i = 1; i <= numv; ++i) {
    j = sta[i];
    sta[i] = toc;
    toc += j;
  }

  vte = new size_t[toc];

  // Step 3:  Finish constructing the vertex to element list.

  for (i = 0; i < nume; ++i) {
    vtc = elementArray[i].vertex_count();
    vtx = elementArray[i].get_modifiable_vertex_index_array();

    for (j = 0; j < vtc; ++j) {
      vte[sta[vtx[j]]++] = i;
    }
  }

  for (i = numv; i > 0; --i) {
    sta[i] = sta[i-1];
  }
  sta[i] = 0;

  // Step 4:  Begin the reodering by computing the vertex furthest from the
  //          origin.

  max = -1.0;
  idx =  0;

  for (i = 0; i < numv; ++i) {
    val = vertexArray[i].length_squared();
    if (val > max) {
      max = val;
      idx = i+1;
    }
  }

  // Step 5:  Perform a breadth first search to find the ordering.

  memset(per, 0, numv*sizeof(size_t));

  loc = 0;
  while (idx > 0) {
    // The vertex referenced by idx has not been visited yet.  Insert it
    // into the queue for processing.
    --idx;

    q1l = 1;
    que1[0] = idx;
    per[idx] = 1;

    while (q1l) {
      q = 0;
      q2l = 0;

      while (q < q1l) {
        idx = que1[q++];
	ord[loc++] = idx;

        st = sta[idx];
        en = sta[idx+1];
        while (st < en) {
          vtc = elementArray[vte[st]].vertex_count();
          vtx = elementArray[vte[st]].get_modifiable_vertex_index_array();
	  ++st;

          for (j = 0; j < vtc; ++j) {
            idx = vtx[j];
            if (!per[idx]) {
              que2[q2l++] = idx;
              per[idx] = 1;
            }
          }
        }
      }

      q1l = q2l;

      tmp  = que1;
      que1 = que2;
      que2 = tmp;
    }

    if (loc >= numv) {
      break;
    }

    // The mesh is not connected.  There is another piece with some vertices
    // remaining.  Repeat the breadth-first search algorithm on the new
    // mesh component.

    max = -1.0;
    idx =  0;
    for (i = 0; i < numv; ++i) {
      if (!per[i]) {
        val = vertexArray[i].length_squared();
        if (val > max) {
          max = val;
          idx = i+1;
        }
      }
    }
  }

  delete[] que1;
  delete[] que2;

  // Step 6:  Compute the permutation vectors

  pel = new size_t[nume];
  for (i = 0; i < nume; ++i) {
    pel[i] = nume;
  }

  toc = 0;
  for (i = 0; i < numv; ++i) {
    loc = ord[numv-1-i];

    per[loc] = i;

    st = sta[loc];
    en = sta[loc+1];
    while (st < en) {
      loc = vte[st++];
      if (nume == pel[loc]) {
        pel[loc] = toc++;
      }
    }
  }

  delete[] ord;
  delete[] vte;
  delete[] sta;

  // Step 7:  Permute the vertices

  v2a = new MsqVertex[vertexArraySize];
  v2h = new Mesh::VertexHandle[vertexArraySize];

  for (i = 0; i < numv; ++i) {
    loc = per[i];
    v2a[loc] = vertexArray[i];
    v2h[loc] = vertexHandlesArray[i];
  }

  delete[] vertexArray;
  delete[] vertexHandlesArray;

  vertexArray = v2a;
  vertexHandlesArray = v2h;

  // Step 8: Permute the elements and vertex indices for the elements

  e2a = new MsqMeshEntity[elemArraySize];
  e2h = new Mesh::ElementHandle[elemArraySize];

  for (i = 0; i < nume; ++i) {
    vtc = elementArray[i].vertex_count();
    vtx = elementArray[i].get_modifiable_vertex_index_array();

    for (j = 0; j < vtc; ++j) {
      vtx[j] = per[vtx[j]];
    }

    loc = pel[i];
    e2a[loc] = elementArray[i];
    e2h[loc] = elementHandlesArray[i];
  }

  delete[] elementArray;
  delete[] elementHandlesArray;

  elementArray = e2a;
  elementHandlesArray = e2h;

  // Step 9: Finish by deleting allocated memory

  delete[] per;
  delete[] pel;

  // Step 10: Recompute vertex to element mapping if it existed.
 
  if (v2eValid) {
    v2eValid = false;
    generate_vertex_to_element_data();
  }
  return;
}

/*! \fn PatchData::num_free_vertices()
   This function has to iterate through all the PatchData vertices to determine
   the number of free vertices. Use with care ! */
int PatchData::num_free_vertices(MsqError &err)
{
  int num_free_vertices=0;
  
  for (size_t i = 0; i < numVertices; ++i)
    if (vertexArray[i].is_free_vertex())
      ++num_free_vertices;
  
  return num_free_vertices;
}


// #undef __FUNC__
// #define __FUNC__ "PatchData::add_element"
// /*! \fn PatchData::add_element(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh, int* vertex_indices, EntityTopology topo,  MsqError &err)

// \param int* vertex_indices ... those indices corresponds to the indices of
// the element's vertices in the PatchData arrays -- see output
// of the add_vertex function.
// */
// int PatchData::add_element(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh,
//                            size_t* vertex_indices, EntityTopology topo,
//                            MsqError &err)
// {
//   int num_verts = MsqMeshEntity::vertex_count(topo);
//   if (!num_verts)
//     err.set_msg("Attempting to add unknown element type to PatchData.");
//   else if (numElements >= elemArraySize)
//     err.set_msg("No space available. Use reserve_element_capacity().");
//   else
//   {
//       // Set the element's type
//     elementArray[numElements].set_element_type(topo);
//     elementHandlesArray[numElements].mesh = mh;
//     elementHandlesArray[numElements].entity = eh;
//       // Go through each vertex
//     for (int n=0; n<num_verts; ++n)
//     {
//         // Make sure it's a valid index
//       if (vertex_indices[n]>=numVertices)
//         err.set_msg("invalid vertex indices");
//         // Set the element's vertex indices
//       elementArray[numElements].set_vertex_index(n, vertex_indices[n]);
//     }
//     return numElements++;
//   }
//   return -1;
// }

// #undef __FUNC__
// #define __FUNC__ "PatchData::add_triangle"
// /*! \fn PatchData::add_triangle(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh, size_t index_vtx1, size_t index_vtx2, size_t index_vtx3, MsqError &err)

// \brief adds a triangle element to the PatchData object.

// \param int index_vertex1 ... those 3 indices corresponds to the indices of
// the triangle's vertices in the PatchData arrays -- see output
// of the add_vertex function.
// */
// void PatchData::add_triangle(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh,
//                              size_t index_vtx1,
//                              size_t index_vtx2,
//                              size_t index_vtx3,
//                              MsqError &err)
// {
//     // make sure we've got space to add this element
//   if (elemArraySize == numElements)
//   {
//     err.set_msg("No more space in PatchData element array");
//     return;
//   }
  
//     // checks the indices are valid
//   if (index_vtx1>=numVertices || // index_vtx1<0 ||
//       index_vtx2>=numVertices || // index_vtx2<0 ||
//       index_vtx3>=numVertices // || index_vtx3<0
//       )
//     err.set_msg("invalid vertex indices");
  
//   elementHandlesArray[numElements].mesh = mh;
//   elementHandlesArray[numElements].entity = eh;
//   elementArray[numElements].set_element_type(TRIANGLE);
//   elementArray[numElements].set_vertex_index(0, index_vtx1);
//   elementArray[numElements].set_vertex_index(1, index_vtx2);
//   elementArray[numElements].set_vertex_index(2, index_vtx3);
//   ++numElements;
  
//   return;
// }


#undef __FUNC__
#define __FUNC__ "PatchData::move_vertices"
/*! \fn PatchData::move_vertices(Vector3D dk[], int nb_vtx, double step_size, MsqError &err)

   It is often useful to use the create_coords_momento() function before
   calling this function.
   Compile with -DMSQ_DBG1 or higher to check that fixed vertices
   are not moved with that call.

   \param dk: must be a [nb_vtx][3] array of doubles that contains the direction in which
          to move each vertex. Fixed vertices moving direction should be zero.
   \param nb_vtx is the number of vertices to move. must corresponds to the number of
          vertices in the PatchData.
   \param step_size will multiply the moving direction given in dk for each vertex.
  */
void PatchData::move_vertices(Vector3D dk[],
                              size_t nb_vtx,
                              double scale,
                              MsqError &err)
{
  if (nb_vtx != numVertices)
  {
    err.set_msg("argument nb_vtx is not equal to the total nb of vertices.");
    MSQ_CHKERR(err);
    return;
  }

    // This will move ALL vertices, only checking for fixed flags after
    // we are done...shouldn't we check for the fixed flag before moving it?
  size_t m=0;
  for (m=0; m<numVertices; ++m)
  {
    vertexArray[m] += (scale * dk[m]);
  }
  
    // Checks that moving direction is zero for fixed vertices.
  MSQ_DEBUG_ACTION(1,{ for (m=0; m<numVertices; ++m)
  {
    Vector3D zero_3d(0.,0.,0.);
    if (   ! vertexArray[m].is_free_vertex()
           && dk[m] != zero_3d                  )
      err.set_msg("moving a fixed vertex.");
    MSQ_CHKERR(err); }
  });
}


#undef __FUNC__
#define __FUNC__ "PatchData::move_free_vertices_constrained"
/*! \fn PatchData::move_free_vertices_constrained(Vector3D dk[], int nb_vtx, double step_size, MsqError &err)
   PatchData::move_free_vertices_constrained() is similar to the function
   PatchData::move_vertices().  The main difference between the two
   being that for this funciton after the free vertices are moved
   as specified by the search direction (dk) and scale factor (step_size),
   they are then projected onto the appropriate geometry.  The second
   difference being that for this function fixed vertices are not moved
   regardless of their corresponding dk direction.
   It is often useful to use the create_coords_momento() function before
   calling this function.
   Compile with -DMSQ_DBG1 or higher to check that fixed vertices
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
  if (nb_vtx != numVertices)
  {
    err.set_msg("The directional vector must be of length numVertices.");
    MSQ_CHKERR(err);
    return;
  }
  
  MsqFreeVertexIndexIterator free_iter(this, err);
  free_iter.reset();
  while (free_iter.next())
  {
    vertexArray[free_iter.value()] += (step_size * dk[free_iter.value()]);
    snap_vertex_to_domain(free_iter.value(), err);
    MSQ_CHKERR(err);
  }
  
    // Checks that moving direction is zero for fixed vertices.
  MSQ_DEBUG_ACTION(1,{ for (size_t m=0; m<numVertices; ++m) {
    Vector3D zero_3d(0.,0.,0.);
    if (   ! vertexArray[m].is_free_vertex()
           && ( dk[m] != zero_3d && dk[m] != -zero_3d)  ) 
    {
      cout << "dk["<<m<<"]: " << dk[m] << endl;
      err.set_msg("moving a fixed vertex.");
      MSQ_CHKERR(err);} 
  }     
  });
}


#undef __FUNC__
#define __FUNC__ "PatchData::set_free_vertices_constrained"
/*! set_free_vertices_constrained is similar to 
PatchData::move_vertices_constrained except the original vertex positions
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
    std::cout << "nb_vtx: "<< nb_vtx << "   mem num vtx; "
              <<  memento->numVertices << std::endl;
    err.set_msg("The directional vector must be the same length as"
                "the number of vertices in the patch.");
    MSQ_CHKERR(err);
    return;
  }
  
  size_t m=0;
  MsqFreeVertexIndexIterator free_iter(this, err);
  free_iter.reset();
  while (free_iter.next())
  {
    m=free_iter.value();
    vertexArray[m] = memento->vertices[m] + (step_size * dk[m]);
    snap_vertex_to_domain(m, err);
    MSQ_CHKERR(err);
  }
  
    // Checks that moving direction is zero for fixed vertices.
  MSQ_DEBUG_ACTION(3,{ for (m=0; m<numVertices; ++m)
  {
    Vector3D zero_3d(0.,0.,0.);
    if (   ! vertexArray[m].is_free_vertex()
           && ( dk[m] != zero_3d && dk[m] != -zero_3d)  ) 
    {
      cout << "dk["<<m<<"]: " << dk[m] << endl;
      err.set_msg("moving a fixed vertex.");
      MSQ_CHKERR(err);
    }
  }
  });
}


#undef __FUNC__
#define __FUNC__ "PatchData::get_max_vertex_movement_squared"
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
  int m=0;
  Vector3D temp_vec;
  double temp_dist=0.0;
  double max_dist=0.0;
  MsqFreeVertexIndexIterator free_iter(this, err);
  free_iter.reset();
  while (free_iter.next())
  {
    m=free_iter.value();
    temp_vec=vertexArray[m] - memento->vertices[m];
    temp_dist=temp_vec.length_squared();
    if(temp_dist>max_dist)
    {
      max_dist=temp_dist;
    }
    MSQ_CHKERR(err);
  }
  return max_dist;
}

#undef __FUNC__
#define __FUNC__ "PatchData::set_all_vertices_soft_fixed"
/*!
 */
void PatchData::set_all_vertices_soft_fixed(MsqError &err)
{
  for(size_t i=0;i<numVertices;++i)
    vertexArray[i].set_soft_fixed_flag();
}

#undef __FUNC__
#define __FUNC__ "PatchData::set_free_vertices_soft_fixed"
/*!
 */
void PatchData::set_free_vertices_soft_fixed(MsqError &err)
{
  for(size_t i=0;i<numVertices;++i){
    if(vertexArray[i].is_free_vertex())
      vertexArray[i].set_soft_fixed_flag();
  }
}

#undef __FUNC__
#define __FUNC__ "PatchData::set_all_vertices_soft_free"
/*!
 */
void PatchData::set_all_vertices_soft_free(MsqError &/*err*/)
  {
    for(size_t i=0;i<numVertices;++i)
      vertexArray[i].remove_soft_fixed_flag();
  }
  
/*! \fn PatchData::get_element_vertex_coordinates(size_t elem_index, std::vector<Vector3D> &coords, MsqError &err)

    \param elem_index The element index in the Patch
    \param coords This std::vector will have the coordinates appended to it.
    If necessary, make sure to clear the vector before calling the function.
  */
void PatchData::get_element_vertex_coordinates(
  size_t elem_index,
  std::vector<Vector3D> &coords,
  MsqError& /*err*/)
{
    // Check index
  if (elem_index >= numElements)
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
  std::vector<size_t> &vertex_indices,
  MsqError& /*err*/)
{
    // Check index
  if (elem_index >= numElements)
    return;
  
    // Ask the element for its vertex indices
  elementArray[elem_index].get_vertex_indices(vertex_indices);
}


#undef __FUNC__
#define __FUNC__ "PatchData::get_vertex_element_indices" 
void PatchData::get_vertex_element_indices(size_t vertex_index,
                                           std::vector<size_t> &elem_indices,
                                           MsqError &err) 
{
    // Check index
  if (vertex_index >= numVertices)
  {
    return;
  }
  
    // Make sure we've got the data
  if (!v2eValid || !v2eOffset)
  {
    generate_vertex_to_element_data();
  }
  
    // Find the starting point for this vertex's data
  size_t *pos = v2E + v2eOffset[vertex_index];
  size_t elem_count = *pos;
  elem_indices.reserve(elem_indices.size() + elem_count);
    // Add each element index to the list
  while (elem_count--)
    elem_indices.push_back(*(++pos));
}

/*!
    \brief This function fills a std::vector<size_t> with the indices
    to vertices connected to the given vertex by an edge.  If vert_indices
    is not initially empty, the function will not delete the current
    contents.  Instead, it will append the new indices at the end of
    the vector.

*/
#undef __FUNC__
#define __FUNC__ "PatchData::get_adjacent_vertex_indices" 
void PatchData::get_adjacent_vertex_indices(size_t vertex_index,
                                            std::vector<size_t> &vert_indices,
                                            MsqError &err)
{
  generate_vertex_to_element_data();
  
    //First get elems attached to vertex[vertex_index]
  std::vector<size_t> elem_indices;
  std::vector<size_t> temp_vert_indices;
  std::vector<size_t>::iterator iter;
  size_t cur_vert;
  int found=0;
  get_vertex_element_indices(vertex_index, elem_indices,err);
  MSQ_CHKERR(err);
  MsqMeshEntity* elems=get_element_array(err);MSQ_CHKERR(err);
    //get nodes attached to vertex_index... with some duplication
  while(!elem_indices.empty()){
    elems[elem_indices.back()].get_connected_vertices(vertex_index, temp_vert_indices,err); MSQ_CHKERR(err);;
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
                                                std::vector<size_t> &adj_ents,
                                                MsqError &err)
{
  //reset the vector
  adj_ents.clear();
    //This should probably be removed
  generate_vertex_to_element_data();
    //vertices of this entity (given by ent_ind)
  std::vector<size_t> verts;
    //vector to store elements attached to the vertices in verts
  std::vector<size_t> elem_on_vert[MSQ_MAX_NUM_VERT_PER_ENT];
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
#undef __FUNC__
#define __FUNC__ "PatchData::update_mesh" 
void PatchData::update_mesh(MsqError &/*err*/)
{
  if (!meshSet)
    return;

  meshSet->update_mesh(*this);
}

#undef __FUNC__
#define __FUNC__ "PatchData::generate_vertex_to_element_data" 
void PatchData::generate_vertex_to_element_data(size_t num_vertex_uses)
{
  FUNCTION_TIMER_START(__FUNC__);
  
    // Skip if data already exists
  if (v2eValid && v2E && v2eOffset)
  {
    FUNCTION_TIMER_END();
    return;
  }
  
    // Create an array that will temporarily hold the number of
    // times each vertex is used in an element.
  if (v2eOffsetSize < num_vertices())
  {
    v2eOffset = new size_t[num_vertices()];
    v2eOffsetSize = num_vertices();
  }
  
    // Initialize each use count to zero
  memset(v2eOffset, 0, num_vertices()*sizeof(size_t));
  
    // Find out how many vertex uses we've got...
  size_t elem_num, vert_num;
  const size_t* vert_array;
  if (num_vertex_uses == 0)
  {
    for (elem_num = num_elements(); elem_num--; )
    {
      num_vertex_uses += elementArray[elem_num].vertex_count();
      vert_array = elementArray[elem_num].get_vertex_index_array();
        // Go through the element, add to each vertex's use count
      for (vert_num = elementArray[elem_num].vertex_count(); vert_num--; )
      {
        v2eOffset[vert_array[vert_num]]++;
      }
    }
  }
    // For each vertex, we store the number of elements along
    // with the element indices.  Increase array size accordingly.
  num_vertex_uses += num_vertices();
  
    // Now that we know how many vertex uses, allocate if necessary
  if (v2eSize < num_vertex_uses)
  {
    v2E = new size_t[num_vertex_uses];
    v2eSize = num_vertex_uses;
  }
  
    // Convert the uses counts to array offsets.
  elem_num = 0;
  for (vert_num = 0; vert_num < num_vertices(); vert_num++)
  {
    size_t temp = v2eOffset[vert_num] + 1;
    v2eOffset[vert_num] = elem_num;
    v2E[elem_num] = 0;
    elem_num += temp;
  }
  
    // Finally, store the v2E data
  for (elem_num = 0; elem_num < num_elements(); elem_num++)
  {
    vert_array = elementArray[elem_num].get_vertex_index_array();
    for (vert_num = elementArray[elem_num].vertex_count(); vert_num--; )
    {
        // Point at the start of the vertex's data
      size_t* ptr = v2E + v2eOffset[vert_array[vert_num]];
        // add one to the number of elements already stored for this element
      (*ptr)++;
        // Store the element number
      ptr[*ptr] = elem_num;
    }
  }
  
  v2eValid = true;
  
  FUNCTION_TIMER_END();
}


#undef __FUNC__
#define __FUNC__ "PatchData::get_subpatch" 
void PatchData::get_subpatch(size_t center_vertex_index,
                             PatchData &subpatch,
                             MsqError &err)
{
    // Make sure we're in range
  if (center_vertex_index >= numVertices)
  {
    err.set_msg("Invalid index for center vertex");
    return;
  }

    // Make sure we've got the vertex-to-element connectivity
  generate_vertex_to_element_data();

    // Make sure we've got the subpatchIndexArray
  if (subpatchIndexSize < numVertices)
  {
    if (subpatchIndexSize > 0)
      delete [] subpatchIndexArray;
    subpatchIndexArray = new size_t[numVertices];
    subpatchIndexSize = numVertices;
  }
  
    // Get the number of elements attached to the requested vertex
  size_t num_elems = v2E[ v2eOffset[center_vertex_index] ];
  subpatch.reserve_element_capacity(num_elems, err);
  subpatch.numElements = num_elems;
  
    // Loop through each element, count the vertices
  size_t which_elem, which_vert, total_verts = 0;
  for (which_elem = num_elems; which_elem--; )
  {
    size_t elem_index = v2E[v2eOffset[center_vertex_index] + which_elem];
    MsqMeshEntity& elem = elementArray[elem_index];
    for (which_vert = elem.vertex_count(); which_vert--; )
    {
      size_t vert_index = elem.get_vertex_index(which_vert);
      if (subpatchIndexArray[vert_index] == 0)
      {
          // This will give the vert a 1-based index.
        subpatchIndexArray[vert_index] = ++total_verts;
      }
    }
  }
  
    // Now that we know how many verts, allocate...
  subpatch.reserve_vertex_capacity(total_verts, err);
  subpatch.numVertices = total_verts;
  
    // Loop through the elements again, placing verts and elems into subpatch
  total_verts = 1;
  for (which_elem = num_elems; which_elem--; )
  {
    size_t elem_index = v2E[v2eOffset[center_vertex_index] + which_elem];
    MsqMeshEntity& elem = elementArray[elem_index];
    for (which_vert = elem.vertex_count(); which_vert--; )
    {
      size_t vert_index = elem.get_vertex_index(which_vert);

        // If we haven't put this vertex into the subpatch yet...
      if (subpatchIndexArray[vert_index] == total_verts)
      {
          // Add this vertex to the subpatch
        subpatch.vertexArray[total_verts-1] = vertexArray[vert_index];
          // Advance total_verts so we know which vertex to add next
        total_verts++;
      }
      
        // Add this vertex to the new element's array
      subpatch.elementArray[which_elem].
        set_vertex_index(which_vert, subpatchIndexArray[vert_index]-1);
    }
  }

    // Loop one final time, setting all vertex indices back to zero
  for (which_elem = num_elems; which_elem--; )
  {
    size_t elem_index = v2E[v2eOffset[center_vertex_index] + which_elem];
    MsqMeshEntity& elem = elementArray[elem_index];
    for (which_vert = elem.vertex_count(); which_vert--; )
    {
      size_t vert_index = elem.get_vertex_index(which_vert);
      subpatchIndexArray[vert_index] = 0;
    }
  }
}

//! Adjust the position of the specified vertex so that it
//! lies on its constraining domain.  The actual domain constraint
//! is managed by the TSTT mesh implementation
void PatchData::snap_vertex_to_domain(size_t vertex_index, MsqError &err)
{
  if (meshSet && meshSet->get_domain_constraint())
  {
    meshSet->get_domain_constraint()->snap_to(vertexHandlesArray[vertex_index],
                                              vertexArray[vertex_index]);
  }
}


#undef __FUNC__
#define __FUNC__ "PatchData::get_domain_normal_at_vertex"
void PatchData::get_domain_normal_at_vertex(size_t vertex_index,
                                   Vector3D &surf_norm,
                                   MsqError &err) const
{
  if (meshSet && meshSet->get_domain_constraint())
  {
    surf_norm = vertexArray[vertex_index];
    meshSet->get_domain_constraint()->normal_at(
      vertexHandlesArray[vertex_index],
      surf_norm);
  }
  else
    err.set_msg("No domain constraint set.");
}


#undef __FUNC__
#define __FUNC__ "PatchData::get_domain_normal_at_element"
void PatchData::get_domain_normal_at_element(size_t elem_index,
                                             Vector3D &surf_norm,
                                             MsqError &err) const
{
  if (meshSet && meshSet->get_domain_constraint())
  {
    elementArray[elem_index].get_centroid(surf_norm, *this, err); MSQ_CHKERR(err);
    meshSet->get_domain_constraint()->normal_at(
      elementHandlesArray[elem_index],
      surf_norm);
  }
  else
    err.set_msg("No domain constraint set.");
}
