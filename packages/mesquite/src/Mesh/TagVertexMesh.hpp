/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TagVertexMesh.hpp
 *  \brief Definition of Mesquite::TagVertexMesh class
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TAG_VERTEX_MESH_HPP
#define MSQ_TAG_VERTEX_MESH_HPP

#include "Mesquite.hpp"
#include "MeshInterface.hpp"

namespace Mesquite {

/**\brief Store alternate vertex coordinates in tags. 
 *
 * This class implements a decorator pattern for the Mesquite::Mesh
 * interface where alternate vertex coordinates are stored in a tag
 * on the original mesh.  The vertex coordinates are the same as that
 * of the decorated Mesh interface until they are set.  Once the coorindates
 * of a vertex are set using an instance of this class, the modified
 * vertex coordinates will be stored in a tag and returned from subsequent
 * queries of the vertex coordinates through this interface.
 *
 * The tag used to store alternate vertex coordinates is created and set
 * for all vertices when the coordinates of the any vertex are changed.
 * The tag type is a vector of three doubles.  
 */
class MESQUITE_EXPORT TagVertexMesh : public Mesh
{
  private:
  
    msq_std::string tagName; //< Name of tag storing vertex coordinates
    TagHandle tagHandle;     //< Handle of tag storing vertex coordinates
    bool haveTagHandle;      //< True if tagHandle is set
    bool cleanUpTag;         //< If true, destroy tag in destructor
    Mesh* realMesh;          //< The actual Mesh this instance decorates
  
      /**\brief common code for constructor, set_mesh, and set_tag_name */
    void initialize( Mesh* mesh, msq_std::string name, bool init, MsqError& );
      /**\brief copy real coordinate values into tag data */
    void copy_all_coordinates( MsqError& err );
      /**\brief if cleanUpTag, delete tag and clear handle */
    void check_remove_tag( MsqError& err );
  
  public:
  
    /** 
     *\param real_mesh  The mesh from which to aquire topology information
     *                  and vertex coordinates, and upon which to store
     *                  tags.
     *\param initialize If the tag for storing alternate vertex coordinates
     *                  exists on the real mesh, those alternate vertex
     *                  coordinates will be retained and used if this
     *                  is false.
     *\param clean_up_tag_data If true, tag storing alternate vertex
     *                  coordinates will be removed when this object
     *                  is destroyed.
     *\param tag_name Name of tag in which to store alternate vertex coordinates.
     */
    TagVertexMesh( MsqError& err,
                   Mesh* real_mesh,
                   bool initialize = true,
                   bool clean_up_tag_data = true,
                   msq_std::string tag_name = "" );
    
      /** Destroy tag data for alternate coordinates if
       *  clean_up_tag_data is true. 
       */
    virtual ~TagVertexMesh();
    
    /**\brief Change the Mesh instance used as the real mesh. 
     *
     * Change the Mesh instance orignially specified in the 
     * constructor.
     * Note: Calling this function changes the handle space for
     *       mesh entities, invalidating any previous handle values, 
     *       iterators, etc. returned by the class instance.
     * Note: If clean_up_tag_data is true, calling this function
     *       will remove any stored alternate vertex coordinates 
     *        from the previous mesh.
     *\param init If the tag for storing alternate vertex coordinates
     *            exists on the real mesh, those alternate vertex
     *            coordinates will be retained and used if this
     *            is false.
     */
    void set_mesh( Mesh* real_mesh, bool init, MsqError& err );
    
    /**\brief Get the real Mesh instance */
    Mesh* get_mesh( ) const { return realMesh; }
    
    /**\brief Set tag cleanup behavior
     *
     * If true, class will remove any tag data storing alternate
     * vertex coordinates from the real mesh when a) the real Mesh
     * instance is changed or b) this object instance is destroted.
     */
    void should_clean_up_tag_data( bool value ) { cleanUpTag = value; }
    
    /**\brief Will tag storing alternate coordinates be destroyed. */
    bool will_clean_up_tag_data() const { return cleanUpTag; }
    
    /**\brief Get name of tag used to store alternate vertex coordinates. */
    msq_std::string get_tag_name() const { return tagName; }

    /**\brief Set tag name used to store alternate vertex coordinates
     *
     * Change the tag name used to store alternate vertex coordinates.
     * Note:  Changing the tag name will result in the loss of any
     *        alternate vertex coordinates saved using the previous
     *        tag name.
     * Note:  If clean_up_tag_data is true, calling this function
     *        will result in the removal of the previous tag and
     *        any coordinates stored using that tag.
     *\param init  If the new tag already exists, any 
     *        coordinates stored in that tag will be used if this
     *        argument is false.  If this argument is true, the
     *        alternate coordinate values will be initialized to
     *        the true coordinate values in the real Mesh.
     */
    void set_tag_name( msq_std::string name, bool init, MsqError& err );
  
    /**\brief clear all alternate vertex coordinate values
     *
     * Clear all alternate vertex coordinate values and 
     * revert to coordinates as stored in real mesh.
     */
    void clear( MsqError& err );
    
    
//************ Operations on entire mesh ****************

    virtual int get_geometric_dimension(MsqError &err);

    virtual void get_all_elements( msq_std::vector<ElementHandle>& elements,
                                   MsqError& err );

    virtual void get_all_vertices( msq_std::vector<VertexHandle>& vertices,
                                   MsqError& err );

    virtual VertexIterator* vertex_iterator(MsqError &err);

    virtual ElementIterator* element_iterator(MsqError &err);

//************ Vertex Properties ********************

    virtual void vertices_get_fixed_flag( const VertexHandle vert_array[], 
                                          bool fixed_flag_array[],
                                          size_t num_vtx, 
                                          MsqError &err );

    virtual void vertices_get_coordinates( const VertexHandle vert_array[],
                                           MsqVertex* coordinates,
                                           size_t num_vtx,
                                           MsqError &err );

    virtual void vertex_set_coordinates( VertexHandle vertex,
                                         const Vector3D &coordinates,
                                         MsqError &err );

    virtual void vertex_set_byte( VertexHandle vertex,
                                  unsigned char byte, 
                                  MsqError &err);

    virtual void vertices_set_byte( const VertexHandle *vert_array,
                                    const unsigned char *byte_array,
                                    size_t array_size, 
                                    MsqError &err );

    virtual void vertex_get_byte( const VertexHandle vertex,
                                  unsigned char *byte, 
                                  MsqError &err );

    virtual void vertices_get_byte( const VertexHandle *vertex,
                                    unsigned char *byte_array,
                                    size_t array_size, 
                                    MsqError &err );
    
//**************** Vertex Topology *****************    

    virtual void vertices_get_attached_elements( 
                         const VertexHandle* vertex_array,
                         size_t num_vertex,
                         msq_std::vector<ElementHandle>& elements,
                         msq_std::vector<size_t>& offsets,
                         MsqError& err );
    
//*************** Element Topology *************

    virtual void elements_get_attached_vertices(
                                   const ElementHandle *elem_handles,
                                   size_t num_elems,
                                   msq_std::vector<VertexHandle>& vert_handles,
                                   msq_std::vector<size_t>& offsets, 
                                   MsqError &err);
    

    virtual void elements_get_topologies(const ElementHandle *element_handle_array,
                                         EntityTopology *element_topologies,
                                         size_t num_elements, MsqError &err);

//***************  Tags  ***********

    virtual TagHandle tag_create( const msq_std::string& tag_name,
                                  TagType type, unsigned length,
                                  const void* default_value,
                                  MsqError &err);

    virtual void tag_delete( TagHandle handle, MsqError& err );

    virtual TagHandle tag_get( const msq_std::string& name, 
                               MsqError& err );

    virtual void tag_properties( TagHandle handle,
                                 msq_std::string& name_out,
                                 TagType& type_out,
                                 unsigned& length_out,
                                 MsqError& err );

    virtual void tag_set_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       const void* tag_data,
                                       MsqError& err );

    virtual void tag_set_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       const void* tag_data,
                                       MsqError& err );

    virtual void tag_get_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       void* tag_data,
                                       MsqError& err );

    virtual void tag_get_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       void* tag_data,
                                       MsqError& err );

    
//**************** Memory Management ****************

    virtual void release_entity_handles( const EntityHandle *handle_array,
                                         size_t num_handles, 
                                         MsqError &err);

    virtual void release();

};


} // namespace Mesquite

#endif
