// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 16-May-02 at 10:26:21
//  LAST-MOD: 14-Jan-03 at 16:34:34 by Thomas Leurent
//
/*! \file MeshSet.cpp

\brief This files implements all the memory management issues related
to the copy of the original TSTT (or other maybe) mesh entity handles
into Mesquite.
That copy is of course encapsulated in the MeshSet class.
  
    \author Thomas Leurent
    \date 2002-05-16  
 */

#include "TSTT_Base.h"

#include "stdio.h"

#include "MeshSet.hpp"
#include "QualityImprover.hpp"

using namespace Mesquite;  


MeshSet::MeshSet() :
  spaceDim(0),
  fixedVertexTagName("fixed")
{
  currentMesh = meshSet.begin();
  currentVertex = verticesSet.begin();
}

// I have a feeling the destructor will do more later, so
// I didn't put it inline
MeshSet::~MeshSet()
{}


/*! \fn  MeshSet::add_mesh(TSTT::Mesh_Handle mh, MsqError &err)

    adds a TSTT mesh to the MeshSet. If used several times,
    it concatenates the vertices from several TSTT meshes into
    one Mesquite::MeshSet.
  */
#undef __FUNC__
#define __FUNC__ "MeshSet::add_mesh"
void MeshSet::add_mesh(TSTT::Mesh_Handle mh, MsqError &err)
{
  // sets MeshSet::SpaceDim
  int dim;
  TSTT::MeshError info=0;
  TSTT::Mesh_GetGeometricDimension(mh, &dim, &info);
  if (spaceDim == 0) 
    spaceDim = dim;
  else if (dim != spaceDim)
  {
    err.set_msg("Meshes of different dimensions added to the same MeshSet.");
    return;
  }
  
  // adds the TSTT::Mesh_Handle to the MeshSet.
  meshSet.push_front(mh);
}


/*! \fn MeshSet::reset(MsqError &err)

    Resets the MeshSet object.
    The current vertex is set back to the first vertex in the first mesh handle.
  */
#undef __FUNC__
#define __FUNC__ "MeshSet::reset"
void MeshSet::reset(MsqError& /*err*/)
{
  verticesSet.clear();
  currentVertex = verticesSet.begin();
  currentMesh = meshSet.begin();
}


/*! \fn MeshSet::get_next_vertices_set(MsqError &err)

    This function is usually called by MeshSet::get_next_patch().
    It retrieves the next vertices set from the TSTT interface.
    If the end has been reached, it resets the MeshSet and returns false .

    This function should be where memory management is optimized,
    by avoiding to fully duplicate the TSTT information.  

 */
#undef __FUNC__
#define __FUNC__ "MeshSet::get_next_vertices_set"
bool MeshSet::get_next_vertices_set(MsqError &err)
{
  // finds last mesh handle.
  std::list<TSTT::Mesh_Handle>::iterator last_mesh;
  last_mesh = meshSet.end();
  --last_mesh; // we want this to be the last Mesh handle, not one past the end.

  if (meshSet.empty()) {
    err.set_msg("MeshSet does not contain any mesh.");
    return false;
  }
 
  // if MeshSet is in initial state
  if ( verticesSet.empty() )
    currentMesh = meshSet.begin();
  // If we're at the end of the MeshSet::verticesSet but there is more
  else if ( currentVertex==verticesSet.end() && currentMesh!=last_mesh ) {
      currentMesh++;
  }
  else {
    // If we're at the end of the MeshSet
    if ( currentVertex==verticesSet.end() && currentMesh==last_mesh ) {
      // reset MeshSet
      currentMesh = meshSet.begin();
      verticesSet.clear();
      currentVertex = verticesSet.begin();
      return false; 
    }
  }
  
  // we need to create the next vertices set
  int num_vertices=0;
  TSTT::MeshError tstt_err=0;
  TSTT::Entity_Handle *vtx;
  TSTT::Mesh_GetEntities(*currentMesh, TSTT::VERTEX, &vtx, &num_vertices, &tstt_err);
  assert(!tstt_err);        
  /* Copy the vertices handles of one MeshHandle in the MeshSet::verticesSet.
     This can be modified to copy less entities
     and run with smaller merory footprint. */
  for (int i=0; i<num_vertices; ++i)
    {
      verticesSet.push_back(EntityEntry(*currentMesh, vtx[i]));
    }
  currentVertex = verticesSet.begin();
  // Can we free these here?  Are the handles useable in later
  // functions if we free them here?
  TSTT::Mesh_FreeEntityHandles(*currentMesh, vtx, &tstt_err);
  assert(!tstt_err);

  return true;
}


/*! \fn MeshSet::get_next_patch(PatchData &pd, MsqError &err )
  \brief This function fills up a PatchData object with the mesh information
  necessary for optimization algorythms.
  The return value is true as long there exist a next patch, false otherwise.
  The list culling is performed in this function.

  This function is a friend of the PatchData class. Therefore the PatchData
  object passed as an argument is filled up directly.

  \param PatchData &pd: this is the PatchData object that will be filled up.
  \param int num_layers: number of layers of adjacencies included in the
  PatchData object. Set to -1 to retrieve the whole mesh.
*/
#undef __FUNC__
#define __FUNC__ "MeshSet::get_next_patch" 
bool MeshSet::get_next_patch(PatchData &pd, PatchDataParameters &pd_params, MsqError &err )
{

  TSTT::MeshError tstt_err=0;
  PatchData::PatchType patch_type = pd_params.get_patch_type();
  long unsigned int culling_method_bits = pd_params.get_culling_method_bits();

  // *************************************************
  // get the next vertices set from the TSTT interface.
  // This could also be replaced by a TSTT iterator.
  // *************************************************

  
  bool more_vtx;
  // if MeshSet object in initial state
  if (patch_type == PatchData::ELEMENTS_ON_VERTEX_PATCH)
  {
    if ( verticesSet.empty() )
    {
      more_vtx = get_next_vertices_set(err); MSQ_CHKERR(err);
      if ( more_vtx == false )
      {
        return false; 
      }
    }
    else
    {
      currentVertex++;
      // If this is the end of the current vertices list
      if (currentVertex == verticesSet.end() )
      {
        more_vtx = get_next_vertices_set(err);  MSQ_CHKERR(err);
        if ( more_vtx == false )
          return false; // no error needed, we're just at the end of the list.
      }
    }
  }
 
  MSQ_DEBUG_ACTION(3,{std::cout << "  o Patch Type: " << patch_type << std::endl; });
  
  // ***********************************************
  // Gathers the local patch information
  // ***********************************************

  if ( patch_type==PatchData::ELEMENTS_ON_VERTEX_PATCH )
  {
      // checks second argument.
    int num_layers = pd_params.get_nb_layers(err); MSQ_CHKERR(err);
    if (num_layers != 1)
    {
      err.set_msg("no implementation for patch depth !=1 yet."); 
      return false;
    }
    
    TSTT::Entity_Handle* vertices =
       new TSTT::Entity_Handle[MSQ_MAX_NUM_VERT_PER_ENT];
    TSTT::Entity_Handle *patch_elements;
    TSTT::EntityTopology *element_topologies;
    int* csr_pointer; // = new int[500]; //dbg
    int* csr_data; // = new int[500]; //dbg
      //  int num_vertices=0;
    int num_elements;
    Vector3D vertex_coordV;
    
      //   ::set_new_handler(out_of_store);
    
      // This loop is the list culling process
      // also checks if this is the end of the vertices list
      // The culling assumes that you've got the same Mesh_Handle
      // for all vertices in verticesSet, so that the culling
      // tag_handle stays the same.
    bool cull_vertex = true;
    
      // retrieves value of fixedVertexTagName .
    void* bnd_tag_handle;
    if (culling_method_bits & PatchData::NO_BOUNDARY_VTX)
    {
      const char* bnd_tag_name = fixedVertexTagName.c_str();
      TSTT::Mesh_tagGetHandle (currentVertex->mesh, bnd_tag_name,
                               &bnd_tag_handle, &tstt_err);
    }
    while (culling_method_bits!=0 && cull_vertex)
    {
      cull_vertex = false; // don't cull , by default
      
        // Culling of boundary vertices
      if (culling_method_bits & PatchData::NO_BOUNDARY_VTX)
      {
        int* on_boundary = NULL;
        int tag_size = sizeof(int);
        TSTT::Mesh_GetTag_Entity(currentVertex->mesh,
                                 (TSTT::cEntity_Handle)(currentVertex->entity),
                                 bnd_tag_handle, (void**)&on_boundary,
                                 &tag_size, &tstt_err);
        
          // Make sure the call succeeded.
          // NOTE: If the tag doesn't exist, we'll get an error...
        assert(!tstt_err);
        assert(on_boundary != NULL);
        
        if (*on_boundary ==1)
        {
          MSQ_DEBUG_ACTION(2,
          {
            std::cout << "      o Culling vertex "
                      << distance(verticesSet.begin(), currentVertex)
                      << " according to NO_BOUNDARY_VTX." << std::endl; });
          cull_vertex = true;
          currentVertex++;
            // If this is the end of the MeshSet vertices list
          if (currentVertex == verticesSet.end() )
          {
            more_vtx = get_next_vertices_set(err);  MSQ_CHKERR(err);
            if ( more_vtx == false )
              return false; // no error needed, we're just at the end of the list.
          }
        }
          // Delete memory allocated by Mesh_GetTag_Entity.
        delete on_boundary;
      }
      // other type of culling ...
      //else if () {
        //}
    } // end of culling loop
    
      // retrieves array of pointers to regions adjacent to the free vertex...
    num_elements = 0; // this will make TSTT allocate the patch_regions array
    TSTT::Entity_GetAdjacencies(currentVertex->mesh,
                                (TSTT::cEntity_Handle*)&(currentVertex->entity),
                                1, TSTT::REGION, &patch_elements,
                                &csr_pointer, &csr_data, 
                                &num_elements, &tstt_err );
    
      // ... if there are no regions adjacent to the vtx,
      // retrieves adjacent faces.
    if ( num_elements == 0 )
    {
      TSTT::Entity_GetAdjacencies(currentVertex->mesh,
                                  (TSTT::cEntity_Handle*)&(currentVertex->entity),
                                  1, TSTT::FACE, &patch_elements,
                                  &csr_pointer, &csr_data,
                                  &num_elements, &tstt_err);
    }
    assert(!tstt_err);
      //std::cout << "num_elements: " << num_elements << std::endl;
      // If there are neither regions or faces adjacent to the vertex, EXIT_FAILURE.
    if ( num_elements == 0 )
    {
      err.set_msg("no regions or faces are adjacent to free vertex.");
      return false;
    }
    
      // Now, get the types of the adjacent elements.
    element_topologies = new enum TSTT::EntityTopology[num_elements];
    TSTT::Entity_GetTopology(currentVertex->mesh,
                             (TSTT::cEntity_Handle *) patch_elements,
                             num_elements,
                             element_topologies, &tstt_err);
    
      // finds an upper bound to the number of vertices in the patch.
      // The worst case (2 triangles in the patch with center vertex
      // on the boundary) must work.
    int max_num_vertices = 1; // center vtx
    for (int i=0; i<num_elements; ++i)
    {
      int n = MsqMeshEntity::vertex_count(element_topologies[i], err); MSQ_CHKERR(err);
      max_num_vertices += n-1; // no need to count center vertex several times.
    }
    
      // gets rid of previous Patch information (but keeps memory allocated).
    pd.clear();
      // allocates memory for coordinates array
      // if current memory is inapropriate
    pd.reserve_vertex_capacity(max_num_vertices, err); MSQ_CHKERR(err);
    
      // allocates memory for connectivity array
      // if current memory is inapropriate
    pd.reserve_element_capacity(num_elements, err); MSQ_CHKERR(err);
    
    int num_coords=3*MSQ_MAX_NUM_VERT_PER_ENT;
    double* vertex_coords = new double[num_coords];
      // enters the coordinates of the central vertex in PatchData::coordsArray
      // at position 0.
    TSTT::Entity_GetVertexCoords(currentVertex->mesh,
                                 (TSTT::cEntity_Handle *) &(currentVertex->entity),
                                 1, TSTT::INTERLEAVED,
                                 &vertex_coords, &num_coords, &tstt_err);
    pd.add_vertex(currentVertex->mesh, currentVertex->entity,
                  &vertex_coords[0], false, err, MsqVertex::MSQ_NO_VTX_FLAG); MSQ_CHKERR(err);
  
    // For each element within the patch,
    // fill up the adequate info in PatchData.
    for (int e=0; e<num_elements; ++e)
    {
        // gets the vertices of an element ...
      int num_element_vtx=MSQ_MAX_NUM_VERT_PER_ENT;
      TSTT::Entity_GetAdjacencies( currentVertex->mesh,
                                   (TSTT::cEntity_Handle *) &(patch_elements[e]),
                                   1, TSTT::VERTEX, &vertices,
                                   &csr_pointer, &csr_data, 
                                   &num_element_vtx, &tstt_err );
      TSTT::Entity_GetVertexCoords(currentVertex->mesh,
                                   (TSTT::cEntity_Handle *) vertices,
                                   num_element_vtx, TSTT::INTERLEAVED,
                                   &vertex_coords, &num_coords,
                                   &tstt_err);
      assert(!tstt_err);
        // ... enters the coordinates of those vertices in PatchData ...
      size_t vtx_ind[MSQ_MAX_NUM_VERT_PER_ENT];
      for (int n=0; n<num_element_vtx; ++n)
      {
        vtx_ind[n] = pd.add_vertex(currentVertex->mesh, vertices[n],
                                   &vertex_coords[3*n], true, err,
                                   MsqVertex::MSQ_HARD_FIXED); MSQ_CHKERR(err);
      }
        // ... and adds the element to PatchData
      pd.add_element(currentVertex->mesh, patch_elements[e], vtx_ind,
                     tstt_to_mesquite(element_topologies[e]), err); MSQ_CHKERR(err);
    }
    delete[] element_topologies;
    delete[] vertices;
    delete[] vertex_coords;
    
      //cout << "numVertices for local patch: " << pd.num_vertices() << endl;
    
      // TODO : provide patch for several layers of adjacencies.
  
      // free entities allocated in TSTT::Entity_GetAdjacencies()
    
    TSTT::Mesh_FreeEntityHandles(currentVertex->mesh, patch_elements, &tstt_err);
    
    return true;
  }
  
    // ***********************************************
    // Gathers the global patch information
    // ***********************************************
  
  else if (patch_type==PatchData::GLOBAL_PATCH) {
  
    TSTT::Entity_Handle* elem_vtx = new TSTT::Entity_Handle[MSQ_MAX_NUM_VERT_PER_ENT];
    int elem_num_coords=3*MSQ_MAX_NUM_VERT_PER_ENT;
    double* elem_vtx_coords = new double[elem_num_coords];
    TSTT::Entity_Handle *current_mesh_elements;
    TSTT::EntityTopology *element_topologies;
    int* csr_pointer; 
    int* csr_data;
    int num_vertices=0;
    int num_regions=0;
    int num_faces=0;
    int num_elements;
    int num_vertices_current_mesh=0;
    int num_regions_current_mesh=0;
    int num_faces_current_mesh=0;
    int num_elements_current_mesh=0;
    MsqVertex::FlagMaskID fixed_flag;
  
 //   ::set_new_handler(out_of_store);

    // loops over all TSTT mesh handles to find total number of entities
    for ( currentMesh = meshSet.begin();
          currentMesh != meshSet.end();
          ++currentMesh ) {

       // finds the number of vertices in the mesh
       num_vertices_current_mesh=0;
       TSTT::Mesh_Services_GetInt( *currentMesh, "number of vertices",
                                   &num_vertices_current_mesh, &tstt_err);
       num_vertices += num_vertices_current_mesh;

       // finds the number of regions in the mesh
       num_regions_current_mesh=0;
       TSTT::Mesh_Services_GetInt( *currentMesh, "number of regions",
                                   &num_regions_current_mesh, &tstt_err);
       num_regions += num_regions_current_mesh;
       // checks that we do not mix surface meshes with volume meshes
       if (num_regions != 0 && num_regions_current_mesh == 0) {
          err.set_msg("some mesh handles have regions, others do not.");
          MSQ_CHKERR(err);
          return false;
       }

       // finds the number of faces in the mesh
       num_faces_current_mesh=0;
       TSTT::Mesh_Services_GetInt( *currentMesh, "number of faces",
                                   &num_faces_current_mesh, &tstt_err);
       num_faces += num_faces_current_mesh;
    
    }

    // gets rid of previous Patch information (but keeps memory allocated).
    pd.clear();
    // allocates memory for coordinates array
    // if current memory is inapropriate
    pd.reserve_vertex_capacity(num_vertices, err); MSQ_CHKERR(err);
    
    // allocates memory for connectivity array
    // if current memory is inapropriate
    if (num_regions != 0) { 
       num_elements = num_regions; }
    else if (num_faces != 0) {
       num_elements = num_faces; }
    else {
       err.set_msg("no regions or faces found in meshes");
       return false; }
    pd.reserve_element_capacity(num_elements, err); MSQ_CHKERR(err);
    

    // loops over all TSTT mesh handles to fill up PatchData
    for ( currentMesh = meshSet.begin();
          currentMesh != meshSet.end();
          ++currentMesh ) {
      
      // retrieves all the elements of the mesh
      num_elements_current_mesh = 0; // this will make TSTT allocate the patch_regions array
      if (num_regions!=0) {
         TSTT::Mesh_GetEntities(*currentMesh,
                                TSTT::REGION, &current_mesh_elements,
                                &num_elements_current_mesh, &tstt_err );
      }
      // ... if there are no regions,
      // retrieves faces.
      if ( num_regions == 0 ) {
         TSTT::Mesh_GetEntities(*currentMesh,
                                TSTT::FACE, &current_mesh_elements,
                                &num_elements_current_mesh, &tstt_err );
      }
      assert(!tstt_err);

      // Now, get the types of the elements.
      element_topologies = new enum TSTT::EntityTopology[num_elements_current_mesh];
      TSTT::Entity_GetTopology( *currentMesh,
                                (TSTT::cEntity_Handle *) current_mesh_elements, 
                                num_elements_current_mesh,
                                element_topologies, &tstt_err);
      
      // retrieves value of fixedVertexTagName.
      // The value can be different for each TSTT mesh handle.
      const char* bnd_tag_name;
      bnd_tag_name = fixedVertexTagName.c_str(); // converts std::string to C-style string.
      void* bnd_tag_handle;
      TSTT::Mesh_tagGetHandle(*currentMesh, bnd_tag_name, &bnd_tag_handle, &tstt_err);
      
      for (int e=0; e<num_elements_current_mesh; ++e) {
         
         // gets the vertices of an element ...
         int num_element_vtx=MSQ_MAX_NUM_VERT_PER_ENT;
         TSTT::Entity_GetAdjacencies( *currentMesh,
                                      (TSTT::cEntity_Handle *) &(current_mesh_elements[e]),
                                      1, TSTT::VERTEX, &elem_vtx,
                                      &csr_pointer, &csr_data, 
                                      &num_element_vtx, &tstt_err );
         TSTT::Entity_GetVertexCoords(*currentMesh,
                                      (TSTT::cEntity_Handle *) elem_vtx,
                                      num_element_vtx, TSTT::INTERLEAVED,
                                      &elem_vtx_coords, &elem_num_coords,
                                      &tstt_err);
         assert(!tstt_err);

         // ... enters the coordinates of those vertices in PatchData ...
         size_t vtx_ind[MSQ_MAX_NUM_VERT_PER_ENT];
         for (int n=0; n<num_element_vtx; ++n) {
            // retrieves the fixed tag value
            int* on_boundary = NULL;
            int tag_size = sizeof(int);
            fixed_flag = MsqVertex::MSQ_NO_VTX_FLAG;
            TSTT::Mesh_GetTag_Entity(*currentMesh,
                                     (TSTT::cEntity_Handle) elem_vtx[n],
                                     bnd_tag_handle, (void**)&on_boundary,
                                     &tag_size, &tstt_err);
            assert(!tstt_err);
            assert(on_boundary != NULL);
            if (*on_boundary==1) {
               fixed_flag = MsqVertex::MSQ_HARD_FIXED;
            }
            // (this retrieves index of previously added vertices). 
            vtx_ind[n] = pd.add_vertex(*currentMesh, elem_vtx[n],
                                       &elem_vtx_coords[3*n], true, err, fixed_flag);
            MSQ_CHKERR(err);
         }

         // ... and adds the element to PatchData.
         pd.add_element(*currentMesh, current_mesh_elements[e],
                        vtx_ind, tstt_to_mesquite(element_topologies[e]), err);
         MSQ_CHKERR(err);
      }
      
      delete[] element_topologies;
            
        // TODO : provide patch for several layers of adjacencies.
      
        // free entities allocated in TSTT::Entity_GetAdjacencies()
      TSTT::Mesh_FreeEntityHandles(*currentMesh, current_mesh_elements, &tstt_err);
    }
    
    delete[] elem_vtx;
    delete[] elem_vtx_coords;
    return true;
  }  
  
  
  err.set_msg("no implementation for specified patch type.");
  return false;
  
}
