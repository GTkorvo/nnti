/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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

    (2007) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file MsqIMesh.hpp
 *  \brief Adaptor for ITAPS iMesh interface
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_I_MESH_HPP
#define MSQ_I_MESH_HPP

#include "MeshInterface.hpp"
#include "iMesh.h"

namespace MESQUITE_NS {

/** The name of the tag (integer) that Mesquite will use
 *  to store internal data
 */
const char* const VERTEX_BYTE_TAG_NAME  = "MesquiteVertexByte";

/** The name of the tag (integer) Mesquite expects to be non-zero
 *  for vertices which are not to be moved by Mesquite
 */
const char* const VERTEX_FIXED_TAG_NAME = "MesquiteVertexFixed";

/** The name of the tag (integer) Mesquite expects to be non-zero
 *  for vertices that are higher-order nodes slaved to their logical
 *  position.
 */
const char* const VERTEX_SLAVED_TAG_NAME = "MesquiteVertexSlaved";

/**\class MsqIMesh
 *\brief Mesquite iMesh Adapter
 *
 * Adpater for interfacing Mesquite with an application that provides
 * the ITAPS iMesh interface for interacting with mesh data.
 */
class MsqIMesh : virtual public Mesquite::Mesh
{
public:
//********* Functions that are NOT inherited ************

  MsqIMesh();
  virtual ~MsqIMesh();

  MsqIMesh( iMesh_Instance imesh, iBase_EntitySetHandle meshset,
	    iBase_EntityType element_dimension,
            MsqError& err,
	    const char* fixed_tag_name = VERTEX_FIXED_TAG_NAME,
	    const char* slaved_tag_name= VERTEX_SLAVED_TAG_NAME );
  
  MsqIMesh( iMesh_Instance imesh, 
            MsqError& err,
	    const char* fixed_tag_name = VERTEX_FIXED_TAG_NAME,
	    const char* slaved_tag_name= VERTEX_SLAVED_TAG_NAME  );
  
  virtual void init_active_mesh( iMesh_Instance mesh, 
                                 MsqError& err,
				 const char* fixed_tag_name = VERTEX_FIXED_TAG_NAME,
				 const char* slaved_tag_name = VERTEX_SLAVED_TAG_NAME );
    
  /** \brief set mesh to be smoothed.
   *
   * Set the mesh which Mesquite is to smooth.  Optionally
   * specify fixed vertices.
   * NOTE: If an active set is not specified, the default
   *       is to use the global set (the ENTIRE mesh.)
   *
   *\param element_set ITAPS entity set handle for set containing
   *                  mesh elements and vertices for which quality 
   *                  is to be improved.
   */
  virtual void set_active_set( iBase_EntitySetHandle meshset, 
                               iBase_EntityType element_dimension,
                               MsqError& err );
  
  virtual iMesh_Instance get_imesh_instance() const;
  virtual iBase_EntitySetHandle get_entity_set() const;

//********* Functions that ARE inherited ************

      /**\brief Get dimension of vertex coordinates (2D vs. 3D). */
    virtual int get_geometric_dimension(Mesquite::MsqError &/*err*/);
    
    /** \brief Get handles for all elemnents */
    virtual void get_all_elements( std::vector<ElementHandle>& elements, 
                                   MsqError& err );
    
    /** \brief Get handles for all vertices */
    virtual void get_all_vertices( std::vector<VertexHandle>& vertices, 
                                   MsqError& err );

      /**\brief Query "fixed" flag for a vertex */
    virtual void vertices_get_fixed_flag( const VertexHandle vert_array[], 
                                          bool fixed_flag_array[],
                                          size_t num_vtx, 
                                          MsqError &err);

    virtual void vertices_get_slaved_flag( const VertexHandle vert_array[], 
                                           bool slaved_flag_array[],
                                           size_t num_vtx, 
                                           MsqError &err );
 
      /**\brief Get vertex coordinates */
    virtual void vertices_get_coordinates( const VertexHandle vert_array[],
                                           MsqVertex* coordinates,
                                           size_t num_vtx, 
                                           MsqError &err);
      /**\brief Set vertex coordinates */
    virtual void vertex_set_coordinates( VertexHandle vertex,
                                         const Vector3D &coordinates, 
                                         MsqError &err);
    
      /**\brief Set vertex mark */
    virtual void vertex_set_byte( VertexHandle vertex,
                                  unsigned char byte, 
                                  MsqError &err);
      /**\brief Set vertex mark */
    virtual void vertices_set_byte( const VertexHandle *vert_array,
                                    const unsigned char *byte_array,
                                    size_t array_size, 
                                    MsqError &err);
    
      /**\brief Get vertex mark */
    virtual void vertex_get_byte( VertexHandle vertex,
                                  unsigned char *byte, 
                                  MsqError &err);
      /**\brief Get vertex mark */
    virtual void vertices_get_byte( const VertexHandle *vert_array,
                                    unsigned char *byte_array,
                                    size_t array_size, 
                                    MsqError &err);
    
      /**\brief Get vertex-to-element adjacencies */
    virtual void vertices_get_attached_elements( const VertexHandle* vertex_array,
                                                 size_t num_vertices,
                                                 std::vector<ElementHandle>& elements,
                                                 std::vector<size_t>& offsets,
                                                 MsqError& err );
    
      /**\brief Get element connectivity */
    virtual void elements_get_attached_vertices( const ElementHandle *elem_handles,
                                                 size_t num_elems,
                                                 std::vector<VertexHandle>& vertices,
                                                 std::vector<size_t>& offsets,
                                                 MsqError& err );
    
  
      /**\brief Return topology type enum for an array of elements */
    virtual void elements_get_topologies( const ElementHandle *element_handle_array,
                                          EntityTopology *element_topologies,
                                          size_t num_elements, 
                                          MsqError &err );
    
//**************** Memory Management ****************
      /**\brief no-op */ 
    virtual void release_entity_handles( const EntityHandle *handle_array,
                                         size_t num_handles, 
                                         MsqError &err );
    
      // Instead of deleting a Mesh when you think you are done,
      // call release().  In simple cases, the implementation could
      // just call the destructor.  More sophisticated implementations
      // may want to keep the Mesh object to live longer than Mesquite
      // is using it.
    virtual void release();

//*************** Tags  ***********

      /** \brief Create a tag
       *
       * Create a user-defined data type that can be attached
       * to any element or vertex in the mesh.  For an opaque or
       * undefined type, use type=BYTE and length=sizeof(..).
       *
       * \param tag_name  A unique name for the data object
       * \param type      The type of the data
       * \param length    Number of values per entity (1->scalar, >1 ->vector)
       * \param default_value Default value to assign to all entities - may be NULL
       * \return - Handle for tag definition 
       */
    virtual TagHandle tag_create( const std::string& tag_name,
                                  TagType type, unsigned length,
                                  const void* default_value,
                                  MsqError &err);
     
      /** \brief Remove a tag and all corresponding data
       *
       * Delete a tag.
       */
    virtual void tag_delete( TagHandle handle, MsqError& err );
    
    
      /** \brief Get handle for existing tag, by name. */
    virtual TagHandle tag_get( const std::string& name, 
                               MsqError& err );
     
      /** \brief Get properites of tag
       *
       * Get data type and number of values per entity for tag.
       * \param handle     Tag to get properties of.
       * \param name_out   Passed back tag name.
       * \param type_out   Passed back tag type.
       * \param length_out Passed back number of values per entity.
       */
    virtual void tag_properties( TagHandle handle,
                                 std::string& name_out,
                                 TagType& type_out,
                                 unsigned& length_out,
                                 MsqError& err );
    
      /** \brief Set tag values on elements
       * 
       * Set the value of a tag for a list of mesh elements.
       * \param handle     The tag 
       * \param num_elems  Length of elem_array
       * \param elem_array Array of elements for which to set the tag value.
       * \param tag_data   Tag data for each element, contiguous in memory.
       *                   This data is expected to be 

       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_set_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       const void* tag_data,
                                       MsqError& err );

      /** \brief Set tag values on vertices
       * 
       * Set the value of a tag for a list of mesh vertices.
       * \param handle     The tag 
       * \param num_elems  Length of node_array
       * \param node_array Array of vertices for which to set the tag value.
       * \param tag_data   Tag data for each element, contiguous in memory.
       *                   This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_set_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       const void* tag_data,
                                       MsqError& err );
    
    
      /** \brief Get tag values on elements
       * 
       * Get the value of a tag for a list of mesh elements.
       * \param handle     The tag 
       * \param num_elems  Length of elem_array
       * \param elem_array Array of elements for which to get the tag value.
       * \param tag_data   Return buffer in which to copy tag data, contiguous 
       *                   in memory.  This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_get_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       void* tag_data,
                                       MsqError& err );
    
      /** \brief Get tag values on vertices.
       * 
       * Get the value of a tag for a list of mesh vertices.
       * \param handle     The tag 
       * \param num_elems  Length of elem_array
       * \param elem_array Array of vertices for which to get the tag value.
       * \param tag_data   Return buffer in which to copy tag data, contiguous 
       *                   in memory.  This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_get_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       void* tag_data,
                                       MsqError& err );
                                       


  protected:
        
    void set_int_tag( void* tag, void* meshset, int value, MsqError& err );

      /** \brief  Call TSTTM::Arr::getEntArrAdj
       *
       * Common code for \ref vertices_get_attached_elements and 
       * \ref elements_get_attached_vertices
       *
       *\param source      Array of handles of source entities to query from
       *\param num_source  The length of \ref source
       *\param target_type The type of entity to query for
       *\param target      The output list of adjacent entities
       *\param offsets     For each entity in \ref source, the offset in 
       *                   \ref target at which the corresponding adjacent
       *                   entities are stored. (output)
       */
    void get_adjacent_entities( const iBase_EntityHandle* source,
                                size_t num_source,
                                iBase_EntityType target_type,
                                std::vector<EntityHandle>& target,
                                std::vector<size_t>& offsets,
                                MsqError& err );

    /** The IMesh instance */
    iMesh_Instance meshInstance;

  private:
      /** \brief Set tag values */
    void tag_set_data ( TagHandle handle,
                        size_t num_elems,
                        const EntityHandle* handle_array,
                        const void* tag_data,
                        MsqError& err );
    
      /** \brief Get tag values */
    void tag_get_data( TagHandle handle,
                       size_t num_elems,
                       const EntityHandle* handle_array,
                       void* tag_data,
                       MsqError& err );
    
    /** Have mesh */
    bool haveMesh;
    /** ITAPS entity set handle for elements to improve */
    //iBase_EntitySetHandle elementSet;
    /** ITAPS entity set handle for nodes to move */
    //iBase_EntitySetHandle nodeSet;
    /** std::set containing elements in elementSet, used
     *  to constrain vertex->element adjaceny queries to
     *  only those elements that are in the input element set.
     */
    //std::vector<iBase_EntityHandle> inputElements;
    
    /** The type of elements contained in the input element set.
     * Should be one of:
     * - iBase_REGION    - volume elements
     * - iBase_FACE      - face/2d elements
     * - iBase_ALL_TYPES - mixed volume and face elements
     */
    iBase_EntityType inputSetType;
    /** The meshset containing the elements to optimize */
    iBase_EntitySetHandle inputSet;
    
    /** Handle for tag used to hold vertex byte */
    iBase_TagHandle byteTag; 
    /** Tag was created in constructor */
    bool createdByteTag;
    /** Handle for tag used to hold vertex-fixed flag */
    iBase_TagHandle fixedTag;
    /** Fixed tag was created in constructor */
    bool createdFixedTag;
    /** Handle for tag used to hold vertex-slaved flag */
    iBase_TagHandle slavedTag;
    /** Handle for the tag used internally to remove duplicates from lists */
//    TagHandle vertexIndexTag;
    /** vertexIndexTag was created in constructor */
//    bool createdVertexIndexTag;
    /** Dimension is queried once during create and cached */
    int geometricDimension;
    /** Map iMesh_EntityTopology to Mesquite::EntityTopology */
    EntityTopology topologyMap[iMesh_ALL_TOPOLOGIES+1];
};

} // namespace Mesquite

#endif
