/*!
  \file   VertexMover.hpp
  \brief  

  The VertexMover Class is the base class for all the smoothing algorythms 

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_VertexMover_hpp 
#define Mesquite_VertexMover_hpp


#include "Mesquite.hpp"
#include "QualityImprover.hpp"
#include "PatchData.hpp"
#include "ObjectiveFunction.hpp"

namespace Mesquite
{

  /*! \class VertexMover
    Base class for all Vertex Movers.
   */  
  class VertexMover : public QualityImprover 
  {
  protected:
    VertexMover();
  public:
    // virtual destructor ensures use of polymorphism during destruction
    virtual ~VertexMover() { };
    
    virtual void loop_over_mesh(MeshSet &ms, MsqError &err);

  protected:

    virtual void initialize(PatchData &pd, MsqError &err) = 0;
    virtual void cleanup() = 0;
    virtual void optimize_vertex_positions(PatchData &pd, 
                                           MsqError &err) = 0; // modifies the PatchData object

    virtual void initialize_mesh_iteration(PatchData &pd, 
                                         MsqError &err) = 0;
    virtual void terminate_mesh_iteration(PatchData &, 
                                         MsqError &err) = 0;

      //!CHECK FEASIBLE IS NOT YET IMPLEMENTED.
    size_t check_feasible(PatchData &pd, MsqError &err);
    
    ObjectiveFunction* objFunc;
  };

  
#undef __FUNC__
#define __FUNC__ "VertexMover::check_feasible"
/*!
  Takes a PatchData object (by reference) and returns whether the
  patch is within the feasible region, 0, or outside the region, 1.
*/
  inline size_t VertexMover::check_feasible(PatchData &pd, MsqError &err)
  {
    MsqMeshEntity* elems=pd.get_element_array(err);
    size_t num_elements=pd.num_elements();
    std::vector<Vector3D> sample_points;
    Vector3D jacobian_vectors[3];
    short num_jacobian_vectors;
    size_t i =0;
    for(i=0;i<num_elements;++i){
      elems[i].get_sample_points(QualityMetric::ELEMENT_VERTICES,sample_points,err);
      std::vector<Vector3D>::iterator iter=sample_points.begin();
      while(iter!=sample_points.end()){
        elems[i].compute_weighted_jacobian(pd, (*iter),
                                           jacobian_vectors,
                                           num_jacobian_vectors, err);
        if(num_jacobian_vectors==2){
            //2-d not yet implemented
        }
        else if(num_jacobian_vectors==3){
          if(jacobian_vectors[0]%(jacobian_vectors[1]*
                                   jacobian_vectors[2])<=0.0){
            return 1;
          }
        }
        ++iter;
      }
    }
    
    return 0;
  }
    

} // namespace
#endif // Mesquite_VertexMover_hpp
