/*!
  \file   MeshInterface.hpp
  \brief  


  \author Darryl Melander
  \author Thomas Leurent
  \date   2003-04-17
*/
#ifndef MESQUITE_INTERFACE_HPP
#define MESQUITE_INTERFACE_HPP

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "Vector3D.hpp"

namespace Mesquite
{
  class EntityIterator;
  typedef EntityIterator VertexIterator;
  typedef EntityIterator ElementIterator;

  inline size_t vertices_in_topology(EntityTopology);

  /*! \class Mesh
     \brief  A Mesquite::Mesh is a collection of mesh elements which are
     composed of mesh vertices.  Intermediate objects are not accessible
     through this interface (where intermediate objects include things
     like the faces of a hex, or an element's edges).
  */
  class Mesh
  {
  public:
//************ Type Definitions **************
      //! Opaque EntityHandle type.
    typedef void* EntityHandle;
    
      // We typedef specific types of EntityHandles just
      // to make it clear what kind of entity is to be
      // returned or used as a function parameter, but the
      // different handle types are not actually distinct.
    typedef EntityHandle VertexHandle;
    typedef EntityHandle ElementHandle;

//************ Operations on entire mesh ****************
      //! Returns whether this mesh lies in a 2D or 3D coordinate system.
    virtual int get_geometric_dimension(MsqError &err) const = 0;
    
      //! Returns the number of entities of the indicated type.
    virtual size_t get_total_vertex_count(MsqError &err) const = 0;
    virtual size_t get_total_element_count(MsqError &err) const = 0;
    
    /*! \brief Fills array with handles to all vertices
                in the mesh.
         \param array_size Must be at least the number of vertices.
                If less than the mesh number of vertices, an error is set
                (and a partial copy is made). */
    virtual void get_all_vertices(VertexHandle *vert_array,
                            size_t array_size, MsqError &err)=0;
    /*! \brief Fills array with handles to all elements
                in the mesh.
         \param array_size Must be at least the number of elements.
                If less than the mesh number of elements, an error is set
                (and a partial copy is made). */
    virtual void get_all_elements(ElementHandle *elem_array,
                            size_t array_size, MsqError &err)=0;
    
      //! Returns a pointer to an iterator that iterates over the
      //! set of all vertices in this mesh.  The calling code should
      //! delete the returned iterator when it is finished with it.
      //! If vertices are added or removed from the Mesh after obtaining
      //! an iterator, the behavior of that iterator is undefined.
    virtual VertexIterator* vertex_iterator(MsqError &err) = 0;
    
      //! Returns a pointer to an iterator that iterates over the
      //! set of all top-level elements in this mesh.  The calling code should
      //! delete the returned iterator when it is finished with it.
      //! If elements are added or removed from the Mesh after obtaining
      //! an iterator, the behavior of that iterator is undefined.
    virtual ElementIterator* element_iterator(MsqError &err) = 0;

//************ Vertex Properties ********************
      //! Returns true or false, indicating whether the vertex
      //! is allowed to be repositioned.  True indicates that the vertex
      //! is fixed and cannot be moved.  Note that this is a read-only
      //! property; this flag can't be modified by users of the
      //! Mesquite::Mesh interface.
    virtual bool vertex_is_fixed(VertexHandle vertex, MsqError &err) = 0;

      //! Returns true or false, indicating whether the vertex
      //! is on the boundary.  Boundary nodes may be treated as
      //! a special case by some algorithms or culling methods.
      //! Note that this is a read-only
      //! property; this flag can't be modified by users of the
      //! Mesquite::Mesh interface.
    virtual bool vertex_is_on_boundary(VertexHandle vertex, MsqError &err) = 0;
    
      //! Get/set location of a vertex
    virtual void vertex_get_coordinates(VertexHandle vertex,
                                        Vector3D &coordinates,
                                        MsqError &err) = 0;
    virtual void vertex_set_coordinates(VertexHandle vertex,
                                        const Vector3D &coordinates,
                                        MsqError &err) = 0;
    
      //! Each vertex has a byte-sized flag that can be used to store
      //! flags.  This byte's value is neither set nor used by the mesh
      //! implementation.  It is intended to be used by Mesquite algorithms.
      //! Until a vertex's byte has been explicitly set, its value is 0.
    virtual void vertex_set_byte (VertexHandle vertex,
                                  unsigned char byte, MsqError &err) = 0;
    virtual void vertices_set_byte (VertexHandle *vert_array,
                                    unsigned char *byte_array,
                                    size_t array_size, MsqError &err) = 0;
    
      //! Retrieve the byte value for the specified vertex or vertices.
      //! The byte value is 0 if it has not yet been set via one of the
      //! *_set_byte() functions.
    virtual void vertex_get_byte(VertexHandle vertex,
                                 unsigned char *byte, MsqError &err) = 0;
    virtual void vertices_get_byte(VertexHandle *vertex,
                                   unsigned char *byte_array,
                                   size_t array_size, MsqError &err) = 0;
    
//**************** Vertex Topology *****************    
      //! Gets the number of elements attached to this vertex.
      //! Useful to determine how large the "elem_array" parameter
      //! of the vertex_get_attached_elements() function must be.
    virtual size_t vertex_get_attached_element_count(VertexHandle vertex,
                                                     MsqError &err) const = 0;
    
      //! Gets the elements attached to this vertex.
    virtual void vertex_get_attached_elements(VertexHandle vertex,
                                              ElementHandle* elem_array,
                                              size_t sizeof_elem_array,
                                              MsqError &err) = 0;
    
//*************** Element Topology *************
    
      //! Gets the number of vertices in this element.
      //! This data can also be found by querying the
      //! element's topology and getting the number
      //! of vertices per element for that topology type.
    virtual size_t element_get_attached_vertex_count(ElementHandle elem,
                                                     MsqError &err) const = 0;

    /*! \brief  
  Returns the vertices that are part of the topological definition of each
  element in the "elem_handles" array.  
  
  When this function is called, the
  following must be true:
   -# "elem_handles" points at an array of "num_elems" element handles.
   -# "vert_handles" points at an array of size "sizeof_vert_handles"
   -# "csr_data" points at an array of size "sizeof_csr_data"
   -# "csr_offsets" points at an array of size "num_elems+1"
      
  When this function returns, adjacency information will be stored
  in csr format:
    -# "vert_handles" stores handles to all vertices found in one
       or more of the elements.  Each vertex appears only
       once in "vert_handles", even if it is in multiple elements.
    -# "sizeof_vert_handles" is set to the number of vertex
       handles placed into "vert_handles".
    -# "sizeof_csr_data" is set to the total number of vertex uses (for
       example, sizeof_csr_data = 6 in the case of 2 TRIANGLES, even if
       the two triangles share some vertices).
    -# "csr_offsets" is filled such that csr_offset[i] indicates the location
       of entity i's first adjacency in "csr_data".  The number of vertices
       in element i is equal to csr_offsets[i+1] - csr_offsets[i].  For this
       reason, csr_offsets[num_elems] is set to the new value of
       "sizeof_csr_data".
    -# "csr_data" stores integer offsets which give the location of
       each adjacency in the "vert_handles" array.

  As an example of how to use this data, you can get the handle of the first
  vertex in element #3 like this:
   \code VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3] ] ] \endcode

  and the second vertex of element #3 like this:
   \code VertexHandle vh = vert_handles[ csr_data[ csr_offsets[3]+1 ] ] \endcode
*/
    virtual void elements_get_attached_vertices(ElementHandle *elem_handles,
                                                size_t num_elems,
                                                VertexHandle *vert_handles,
                                                size_t &sizeof_vert_handles,
                                                size_t *csr_data,
                                                size_t &sizeof_csr_data,
                                                size_t *csr_offsets,
                                                MsqError &err) = 0;
    
      //! Identifies the vertices attached to this element by returning
      //! each vertex's global index.  The vertex's global index indicates
      //! where that element can be found in the array returned by
      //! Mesh::get_all_vertices.
    virtual void element_get_attached_vertex_indices(ElementHandle element,
                                                     size_t *index_array,
                                                     size_t array_size,
                                                     MsqError &err) = 0;
    
      //! Returns the topology of the given entity.
    virtual EntityTopology element_get_topology(ElementHandle entity_handle,
                                              MsqError &err) const = 0;
      //! Returns the topologies of the given entities.  The "entity_topologies"
      //! array must be at least "num_elements" in size.
    virtual void elements_get_topologies(ElementHandle *element_handle_array,
                                         EntityTopology *element_topologies,
                                         size_t num_elements, MsqError &err) = 0;
    
//**************** Memory Management ****************
      //! Tells the mesh that the client is finished with a given
      //! entity handle.  
    virtual void release_entity_handles(EntityHandle *handle_array,
                                        size_t num_handles, MsqError &err) = 0;
    
      //! Instead of deleting a Mesh when you think you are done,
      //! call release().  In simple cases, the implementation could
      //! just call the destructor.  More sophisticated implementations
      //! may want to keep the Mesh object to live longer than Mesquite
      //! is using it.
    virtual void release() = 0;
    
  protected:
      //! Don't allow a Mesh to be deleted directly.
    virtual ~Mesh()
      {}
  };
  
  /*! \class EntityIterator 
  \brief   Iterates through a set of entities.  An EntityIterator 
  is typically obtained via Mesh::vertex_iterator() or
  Mesh::element_iterator().  An iterator obtained in this
  way iterates over the set of all vertices/elements in
  the Mesh from which the iterator was obtained. */
  class EntityIterator
  {
  public:
    virtual ~EntityIterator()
      {}

      //! Moves the iterator back to the first
      //! entity in the list.
    virtual void restart() = 0;
    
      //! *iterator.  Return the handle currently
      //! being pointed at by the iterator.
    virtual Mesh::EntityHandle operator*() const = 0;
    
      //! ++iterator
    virtual void operator++() = 0;
    
      //! Returns false until the iterator has
      //! been advanced PAST the last entity.
      //! Once is_at_end() returns true, *iterator
      //! returns 0.
    virtual bool is_at_end() const = 0;
  };
  

  /*! \class MeshDomain
      The MeshDomain class provides geometrical information concerning the Mesh.
      It is called during surface meshes optimization to figure out the surface normal,
      how to snap vertices back to the surface, etc... . 
    */
  class MeshDomain
  {
  public:
    virtual ~MeshDomain()
      {}

      //! Modifies "coordinate" so that it lies on the
      //! domain to which "entity_handle" is constrained.
      //! The handle determines the domain.  The coordinate
      //! is the proposed new position on that domain.
    virtual void snap_to(Mesh::EntityHandle entity_handle,
                         Vector3D &coordinate) = 0;
    
      //! Returns the normal of the domain to which
      //! "entity_handle" is constrained.  For non-planar surfaces,
      //! the normal is calculated at the point on the domain that
      //! is closest to the passed in value of "coordinate".  If the
      //! domain does not have a normal, or the normal cannot
      //! be determined, "coordinate" is set to (0,0,0).  Otherwise,
      //! "coordinate" is set to the domain's normal at the
      //! appropriate point.
      //! In summary, the handle determines the domain.  The coordinate
      //! determines the point of interest on that domain.
    virtual void normal_at(Mesh::EntityHandle entity_handle,
                           Vector3D &coordinate) = 0;
  };
}

#include "MsqMeshEntity.hpp"

inline size_t Mesquite::vertices_in_topology(Mesquite::EntityTopology topo)
{
  return Mesquite::MsqMeshEntity::vertex_count(topo);
}

#endif
