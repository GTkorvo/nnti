/*!
  \file   PatchData.cpp

  \author Thomas Leurent
  \date   2002-01-17
*/

#include "PatchData.hpp"
#include "MsqMeshEntity.hpp"
#include "MsqFreeVertexIndexIterator.hpp"

#include <list>

using namespace Mesquite;  

#undef __FUNC__
#define __FUNC__ "PatchData::PatchData"
PatchData::PatchData()
  : numVertices(0),
    numElements(0),
    vertexArray(NULL),
    vertexHandlesArray(NULL),
    elementArray(NULL),
    elementHandlesArray(NULL),
    vertexArraySize(0),
    elemArraySize(0),
    elemsInVertex(NULL),
    vertexToElemOffset(NULL)
{}

// // copy function used in copy constructor and assignement
// #undef __FUNC__
// #define __FUNC__ "PatchData::copy" 
// void PatchData::copy(const PatchData &A)
// {
//   storage = A.storage;
//   spaceDim = A.spaceDim;
//   numVertices = A.numVertices;
//   numElements = A.numElements;
//   numFreeVertices = A.numFreeVertices;
//   numFreeElements = A.numFreeElements;

  
//   switch(storage) {
//   case RAW_ARRAYS:
//     // allocates memory and copy the coordinates array
//     coordsArray = new Vector3D[coordsArraySize = A.coordsArraySize];
//     for (int i=0; i<coordsArraySize; i++)
//       coordsArray[i] = A.coordsArray[i];

//     // allocates memory and copy the connectivity array
//     connectivityArray = new ConnectivityArrayT[connectArraySize = A.connectArraySize];
//     for (int i=0; i<connectArraySize; i++)
//       for (int j=0; j<30; j++) {
//         connectivityArray[i].indices[j] = A.connectivityArray[i].indices[j];
//         connectivityArray[i].entity_type = A.connectivityArray[i].entity_type;
//       }

//     break;
//   case MESQUITE_OBJECTS:
//     // allocates memory and copy the vertex pointers array
//     vertexsArray = new MsqVertex*[coordsArraySize = A.coordsArraySize];
//     for (int i=0; i<coordsArraySize; i++)
//       vertexsArray[i] = A.vertexsArray[i];

//     // allocates memory and copy the elements pointers array
//     elementsArray = new MsqMeshEntity*[connectArraySize = A.connectArraySize];
//     for (int i=0; i<connectArraySize; i++)
//       elementsArray[i] = A.elementsArray[i];

//     break;
//   }
// }

// // copy constructor
// #undef __FUNC__
// #define __FUNC__ "PatchData::PatchData(const PatchData )" 
// PatchData::PatchData(const PatchData &A)
// {
// #ifdef MSQ_DEBUG4
//   std::cout << "MSQ_DEBUG4: Executing PatchData copy constructor \n";
// #endif
//   this->copy(A);
// }

// // assignement
// // (performance could be improved by reassigning memory when needed only)
// #undef __FUNC__
// #define __FUNC__ "PatchData::operator=" 
// PatchData& PatchData::operator= (const PatchData &A)
// {
// #ifdef MSQ_DEBUG4
//   std::cout << "MSQ_DEBUG4: Executing PatchData copy assignement \n";
// #endif
//   if (this != &A) {  // beware of self-assignement A=A
//     std::cout<< "deleting coordsarray\n";
//     delete[] coordsArray;
//     std::cout<< "deleting coordsarray\n";
//     delete[] connectivityArray;
//     std::cout<< "deleting connectsarray\n";
//     delete[] vertexsArray;
//     std::cout<< "deleting vertexsarray\n";
//     delete[] elementsArray;
//     std::cout<< "deleting elementarray\n";
//     this->copy(A);
//   }
//   return *this;
// }

// Destructor
#undef __FUNC__
#define __FUNC__ "PatchData::~PatchData" 
PatchData::~PatchData()
{
  delete [] vertexArray;
  delete [] vertexHandlesArray;
  delete [] elementArray;
  delete [] elementHandlesArray;
  delete [] elemsInVertex;
  delete [] vertexToElemOffset;
}


/*! \fn PatchData::num_free_vertices()
   This function has to iterate through all the PatchData vertices to determine
   the number of free vertices. Use with care ! */
int PatchData::num_free_vertices(MsqError &err)
{   
  int num_free_vertices=0;
  MsqFreeVertexIndexIterator ind(this, err);
  ind.reset();
  while (ind.next()) {
    ++num_free_vertices;
  }
  return num_free_vertices;
}


#undef __FUNC__
#define __FUNC__ "PatchData::add_element"
/*! \fn PatchData::add_element(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh, int* vertex_indices, EntityTopology topo,  MsqError &err)

\param int* vertex_indices ... those indices corresponds to the indices of
the element's vertices in the PatchData arrays -- see output
of the add_vertex function.
*/
int PatchData::add_element(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh,
                            int* vertex_indices, EntityTopology topo,
                            MsqError &err)
{
  int num_verts = MsqMeshEntity::vertex_count(topo);
  if (!num_verts)
    err.set_msg("Attempting to add unknown element type to PatchData.");
  else if (numElements >= elemArraySize)
    err.set_msg("No space available. Use reserve_element_capacity().");
  else
  {
      // Set the element's type
    elementArray[numElements].set_element_type(topo);
    elementHandlesArray[numElements].mesh = mh;
    elementHandlesArray[numElements].entity = eh;
      // Go through each vertex
    for (int n=0; n<num_verts; ++n)
    {
        // Make sure it's a valid index
      if ( vertex_indices[n]>=numVertices || vertex_indices[n]<0 )
        err.set_msg("invalid vertex indices");
        // Set the element's vertex indices
      elementArray[numElements].set_vertex_index(n, vertex_indices[n]);
    }
    return numElements++;
  }
  return -1;
}

#undef __FUNC__
#define __FUNC__ "PatchData::add_triangle"
/*! \fn PatchData::add_triangle(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh, size_t index_vtx1, size_t index_vtx2, size_t index_vtx3, MsqError &err)

\brief adds a triangle element to the PatchData object.

\param int index_vertex1 ... those 3 indices corresponds to the indices of
the triangle's vertices in the PatchData arrays -- see output
of the add_vertex function.
*/
void PatchData::add_triangle(TSTT::Mesh_Handle mh, TSTT::Entity_Handle eh,
                             size_t index_vtx1,
                             size_t index_vtx2,
                             size_t index_vtx3,
                             MsqError &err)
{
    // make sure we've got space to add this element
  if (elemArraySize == numElements)
  {
    err.set_msg("No more space in PatchData element array");
    return;
  }
  
    // checks the indices are valid
  if (index_vtx1>=numVertices || // index_vtx1<0 ||
      index_vtx2>=numVertices || // index_vtx2<0 ||
      index_vtx3>=numVertices // || index_vtx3<0
      )
    err.set_msg("invalid vertex indices");
  
  elementHandlesArray[numElements].mesh = mh;
  elementHandlesArray[numElements].entity = eh;
  elementArray[numElements].set_element_type(TRIANGLE);
  elementArray[numElements].set_vertex_index(0, index_vtx1);
  elementArray[numElements].set_vertex_index(1, index_vtx2);
  elementArray[numElements].set_vertex_index(2, index_vtx3);
  ++numElements;
  
  return;
}


#undef __FUNC__
#define __FUNC__ "PatchData::move_free_vertices"
/*! \fn PatchData::move_free_vertices(Vector3D dk[], int nb_vtx, double step_size, MsqError &err)

   It is often useful to use the create_coords_momento() function before
   calling this function.

   \param dk: must be a [nb_vtx][3] array of doubles that contains the direction in which
          to move each free vertex.
   \param nb_vtx is the number of free vertices to move. must corresponds to the number of
          free vertices in the PatchData.
   \param step_size will multiply the moving direction given in dk for each vertex.
  */
void PatchData::move_free_vertices(Vector3D dk[], int nb_vtx,
                                   double step_size, MsqError &err)
{
  if (nb_vtx > numVertices) {
    err.set_msg("argument nb_vtx is higher than the total nb of vertices.");
    MSQ_CHKERR(err);
    return;
  }

  int i=0;
  MsqFreeVertexIndexIterator ind(this, err); MSQ_CHKERR(err);
  ind.reset();
  while (ind.next()) {
      vertexArray[ind.value()] += (step_size * dk[i]);
      ++i;
  }

  if (i != nb_vtx) {
    err.set_msg("Argument nb_vtx is not equal to the number of free vertices.");
    MSQ_CHKERR(err);
  }
}

/*! \fn PatchData::get_element_vertex_coordinates(size_t elem_index, std::vector<Vector3D> &coords, MsqError &err)

    \param elem_index The element index in the Patch
    \param coords This std::vector will have the coordinates appended to it.
    If necessary, make sure to clear the vector before calling the function. 
  */
void PatchData::get_element_vertex_coordinates(
  size_t elem_index,
  std::vector<Vector3D> &coords,
  MsqError &err)
{
    // Check index
  if (elem_index >= numElements)
    return;

    // Ask the element for its vertex indices
  std::vector<size_t> vertex_indices;
  elementArray[elem_index].get_vertex_indices(vertex_indices);
    // Get the coords for each indicated vertex
  for (int i = 0; i < vertex_indices.size(); i++)
    coords.push_back(Vector3D(vertexArray[vertex_indices[i]]));
}

void PatchData::get_element_vertex_indices(
  size_t elem_index,
  std::vector<size_t> &vertex_indices,
  MsqError &err)
{
    // Check index
  if (elem_index >= numElements)
    return;
  
    // Ask the element for its vertex indices
  elementArray[elem_index].get_vertex_indices(vertex_indices);
}

void PatchData::get_vertex_element_indices(
  size_t vertex_index,
  std::vector<size_t> &elem_indices)
{
    // Check index
  if (vertex_index >= numVertices)
    return;
  
    // Make sure we've got the data
  if (!vertexToElemOffset)
    return;
  
    // Find the starting point for this vertex's data
  size_t *pos = elemsInVertex + vertexToElemOffset[vertex_index];
  size_t elem_count = *pos;
  
    // Add each element index to the list
  while (elem_count--)
    elem_indices.push_back(*(++pos));
}


/*! \fn PatchData::update_mesh(MsqError &err)

    \brief This function copies to the TSTT mesh  the changes made to the
    free vertices / elements of the PatchData object.

*/
#undef __FUNC__
#define __FUNC__ "PatchData::update_mesh" 
void PatchData::update_mesh(MsqError &err)
{
  double coordsc[3];
  TSTT::cMesh_Handle mh;
  TSTT::Entity_Handle vertex;
  TSTT::MeshError tstt_err=0;

  for (int n=0; n<numVertices; ++n) {
    mh = vertexHandlesArray[n].mesh;
    vertex = vertexHandlesArray[n].entity;
    
    coordsc[0] = vertexArray[n][0];
    coordsc[1] = vertexArray[n][1];
    coordsc[2] = vertexArray[n][2];
    

    // TODO: creates arrays and pass them in one function call.
    TSTT::Entity_SetVertexCoords( mh, (TSTT::cEntity_Handle*) &vertex,
                                  1, TSTT::INTERLEAVED,
                                  3, coordsc, &tstt_err);
  }
}

void PatchData::generate_vertex_to_element_data()
{
    // Skip if data already exists
  if (elemsInVertex || vertexToElemOffset)
    return;
  
    // Create an array of linked lists, one list per vertex
    // The lists hold element indices.
  std::list<size_t> *element_indices =
    new std::list<size_t>[num_vertices()];
  
    // Go through each element
  size_t total_entries = num_vertices();
  for (int elem_num = num_elements();
       elem_num--; )
  {
    for (int which_node = elementArray[elem_num].vertex_count();
         which_node--; )
    {
        // Put the element number into the vertex's list
      element_indices[elementArray[elem_num].get_vertex_index(which_node)].push_back(elem_num);
      total_entries++;
    }
  }
  
    // Transfer the data in the lists into array data
  vertexToElemOffset = new size_t[num_vertices()];
  elemsInVertex = new size_t[total_entries];
  size_t cur_index = 0;
  for (int i = 0; i < num_vertices(); i++)
  {
      // Save the start point for this vertex's data
    vertexToElemOffset[i] = cur_index;
      // Indicate how many elements in this vertex
    elemsInVertex[cur_index++] = element_indices[i].size();
      // Save the index of each element in this vertex
    for (std::list<size_t>::const_iterator iter = element_indices[i].begin();
         iter != element_indices[i].end();
         ++iter)
    {
      elemsInVertex[cur_index++] = *iter;
    }
  }
  
    // Cleanup
  delete [] element_indices;
}
