/*!
  \file   QualityMetric.cpp
  \brief  

  \author Michael Brewer
  \author Thomas Leurent
  \date   2002-05-14
*/

#include "QualityMetric.hpp"
#include "MsqVertex.hpp"
#include "MsqMeshEntity.hpp"
#include "MsqMessage.hpp"
#include "MsqTimer.hpp"
#include "PatchData.hpp"

using namespace Mesquite;
using std::cout;
using std::endl;
using std::cerr;


#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_element_hessian"
/*!
  \param pd: PatchData that contains the element which Hessian we want.

  \param el: this is the element for which the Hessian will be returned.
  
  \param free_vtces: base address of an array of pointers to the element
  vertices which are considered free for purposes of computing the hessian.
  The vertices within this array must be ordered in the same order
  as the vertices within the element, el. 
  Only the Hessian entries corresponding to a pair of free vertices
  will be non-zero.

  \param grad_vec: this is an array of nve Vector3D, where nve is the total
  number of vertices in the element. Only the entries corresponding to free
  vertices specified in free_vtces will be non-zero. The order is the same
  as the order of the vertices in el.

  \param hessian: this is a 1D array of Matrix3D that will contain the upper
  triangular part of the Hessian. It has size nve*(nve+1)/2, i.e. the number
  of entries in the upper triangular part of a nve*nve matrix.

  \param num_free_vtx: is the number df free vertices in the element.
  Essentially, this gives the size of free_vtces[].  The gradient array has
  the size of the number of vertices in the element, regardless.

  \param metric_value: Since the metric is computed, we return it.
  
  \return true if the element is valid, false otherwise. 
*/
bool QualityMetric::compute_element_hessian(PatchData &pd,
                                            MsqMeshEntity* el,
                                            MsqVertex* free_vtces[],
                                            Vector3D grad_vec[],
                                            Matrix3D hessian[],
                                            int num_free_vtx,
                                            double &metric_value,
                                            MsqError &err)
{
  // first, checks that free vertices order is consistent with the
  // element order. 
  std::vector<size_t> elem_vtx_indices;
  std::vector<size_t>::const_iterator v;
  el->get_vertex_indices(elem_vtx_indices);
  int i;
  v=elem_vtx_indices.begin();
  for (i=0; i<num_free_vtx; ++i) {
    while ( *v != pd.get_vertex_index(free_vtces[i]) ) {
      if ( v==elem_vtx_indices.end() ) {
        err.set_msg("free vertices cannot be given in a different"
                    "order than the element's.");
        return false;
      }
      else  ++v;
    }
  }
    
    
  bool ret=false;
  switch(hessianType)
    {
    case NUMERICAL_HESSIAN:
      ret = compute_element_numerical_hessian(pd, el, free_vtces, grad_vec, hessian,
                                              num_free_vtx, metric_value, err);
      MSQ_CHKERR(err);
      break;
    case ANALYTICAL_HESSIAN:
      ret = compute_element_analytical_hessian(pd, el, free_vtces, grad_vec, hessian,
                                               num_free_vtx, metric_value, err);
      MSQ_CHKERR(err);
      break;
    }
  return ret;
}
   
   
#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_vertex_analytical_gradient"
/*! If that function is not over-riden in the concrete class, the base
    class function makes it default to a numerical gradient.
    \param vertex  Vertex which is considered free for purposes of computing the gradient.
    \param grad_vec Vector where the gradient is stored.
    \param metric_value Since the metric is computed, we return it. 
    \return true if the element is valid, false otherwise.
*/
bool QualityMetric::compute_vertex_analytical_gradient(PatchData &pd,
                                                       MsqVertex &vertex,
                                                       MsqVertex* free_vtces[],
                                                       Vector3D grad_vec[],
                                                       int num_vtx,
                                                       double &metric_value,
                                                       MsqError &err)
{
  PRINT_WARNING("QualityMetric has no analytical gradient defined. ",
                "Defaulting to numerical gradient.\n");
  set_gradient_type(NUMERICAL_GRADIENT);
  return compute_vertex_numerical_gradient(pd, vertex, free_vtces, grad_vec,
                                           num_vtx, metric_value, err);
}

#undef __FUNC__
#define __FUNC__ "QualityMetric::change_metric_type"
void QualityMetric::change_metric_type(MetricType /*t*/, MsqError &err)
{
  err.set_msg("This QualityMetric's MetricType can not be changed.");
}


#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_element_analytical_gradient"
/*! If that function is not over-riden in the concrete class, the base

    Parameters description, see QualityMetric::compute_element_gradient() .

    \return true if the element is valid, false otherwise.
*/
bool QualityMetric::compute_element_analytical_gradient(PatchData &pd,
                                             MsqMeshEntity* element,
                                             MsqVertex* free_vtces[], Vector3D grad_vec[],
                                             int num_free_vtx, double &metric_value,
                                             MsqError &err)
{
  PRINT_WARNING("QualityMetric has no analytical gradient defined. ",
                "Defaulting to numerical gradient.\n");
  set_gradient_type(NUMERICAL_GRADIENT);
  return compute_element_numerical_gradient(pd, element, free_vtces, grad_vec, num_free_vtx, metric_value, err);
}


#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_element_analytical_hessian"
/*! If that function is not over-riden in the concrete class, the base
  class function makes it default to a numerical hessian.
  
  For parameters description, see QualityMetric::compute_element_hessian() .
  
  \return true if the element is valid, false otherwise. 
*/
bool QualityMetric::compute_element_analytical_hessian(PatchData &pd,
                                             MsqMeshEntity* element,
                                             MsqVertex* free_vtces[], Vector3D grad_vec[],
                                             Matrix3D hessian[],
                                             int num_free_vtx, double &metric_value,
                                             MsqError &err)
{
  PRINT_WARNING("QualityMetric has no analytical hessian defined. ",
                "Defaulting to numerical hessian.\n");
  set_hessian_type(NUMERICAL_HESSIAN);
  return compute_element_numerical_hessian(pd, element, free_vtces, grad_vec,
                                           hessian, num_free_vtx, metric_value, err);
}


#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_element_gradient_expanded"
/*!
  Note that for this function, grad_vec should be an array of size the
  number of vertices in el, not of size num_free_vtx.
*/
bool QualityMetric::compute_element_gradient_expanded(PatchData &pd,
                                                      MsqMeshEntity* el,
                                                      MsqVertex* free_vtces[],
                                                      Vector3D grad_vec[],
                                                      int num_free_vtx,
                                                      double &metric_value,
                                                      MsqError &err)
{
  int i, g, e;
  bool ret;
  Vector3D* grad_vec_nz = new Vector3D[num_free_vtx];
  ret = compute_element_gradient(pd, el, free_vtces, grad_vec_nz,
                                 num_free_vtx, metric_value, err);
  MSQ_CHKERR(err);

  std::vector<size_t> gv_i;
  gv_i.reserve(num_free_vtx);
  i=0;
  for (i=0; i<num_free_vtx; ++i) {
    gv_i.push_back( pd.get_vertex_index(free_vtces[i]) );
  }
     
  std::vector<size_t> ev_i;
  el->get_vertex_indices(ev_i);

  bool inc;
  std::vector<size_t>::iterator ev;
  std::vector<size_t>::iterator gv;
  for (ev=ev_i.begin(), e=0; ev!=ev_i.end(); ++ev, ++e) {
    inc = false; g=0;
    gv = gv_i.begin();
    while (gv!=gv_i.end()) {
      if (*ev == *gv) {
        inc = true;
        break;
      }
      ++gv; ++g;
    }
    if (inc == true)
      grad_vec[e] = grad_vec_nz[g];
    else
      grad_vec[e] = 0;
  }
  
  delete []grad_vec_nz;
  return ret;
}
   
   
#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_element_numerical_gradient"
/*!
  Parameters description, see QualityMetric::compute_element_gradient() .
  
  \return true if the element is valid, false otherwise.
*/
bool QualityMetric::compute_element_numerical_gradient(PatchData &pd,
                                             MsqMeshEntity* element,
                                             MsqVertex* free_vtces[],
                                                       Vector3D grad_vec[],
                                             int num_free_vtx, double &metric_value,
                                             MsqError &err)
{
  FUNCTION_TIMER_START(__FUNC__);
    /*!TODO: (MICHAEL)  Try to inline this function (currenlty conflicts
      with MsqVertex.hpp).*/    
  MSQ_DEBUG_PRINT(3,"Computing Numerical Gradient\n");
  
  bool valid=this->evaluate_element(pd, element, metric_value, err); MSQ_CHKERR(err);

  if (!valid)
    return false;
  double delta = 10e-6;
  int counter=0;
  double metric_value1=0;
  for (int v=0; v<num_free_vtx; ++v) 
  {
    /* gradient in the x, y, z direction */
    for (int j=0;j<3;++j) 
    {
        //re-initialize variables.
      valid=false;
      delta = 10e-6;
      counter=0;
        //perturb the node and calculate gradient.  The while loop is a
        //safety net to make sure the epsilon perturbation does not take
        //the element out of the feasible region.
      while(!valid && counter<10){
          // perturb the coordinates of the free vertex in the j direction
          // by delta       
        (*free_vtces[v])[j]+=delta;
          //compute the function at the perturbed point location
        valid=this->evaluate_element(pd, element,  metric_value1, err);
        MSQ_CHKERR(err);
          //compute the numerical gradient
        grad_vec[v][j]=(metric_value1-metric_value)/delta;
          // put the coordinates back where they belong
        (*free_vtces[v])[j] -= delta;
        ++counter;
        delta/=10.0;
      }
      if(counter>=10){
        err.set_msg("Perturbing vertex by delta caused an inverted element.");
      }
      
    }
  }
  FUNCTION_TIMER_END();
  return true;
}


#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_element_numerical_hessian"
/*!
  Note that for this function, grad_vec should be an array of size the
  number of vertices in el, not of size num_free_vtx. Entries that do not correspond
  with the vertices argument array will be null.

  For parameters description, see QualityMetric::compute_element_hessian() .
  
  \return true if the element is valid, false otherwise. 
*/
bool QualityMetric::compute_element_numerical_hessian(PatchData &pd,
                                             MsqMeshEntity* element,
                                             MsqVertex* free_vtces[],
                                             Vector3D grad_vec[],
                                             Matrix3D hessian[],
                                             int num_free_vtx, double &metric_value,
                                             MsqError &err)
{
  FUNCTION_TIMER_START(__FUNC__);
  MSQ_DEBUG_PRINT(3,"Computing Numerical Hessian\n");
  
  bool valid=this->compute_element_gradient_expanded(pd, element, free_vtces, grad_vec,
                                    num_free_vtx, metric_value, err); MSQ_CHKERR(err);
  
  if (!valid)
    return false;
  
  double delta = 10e-6;
  double vj_coord;
  short nve = element->vertex_count();
  Vector3D* grad_vec1 = new Vector3D[nve];
  Vector3D fd;
  std::vector<size_t> ev_i;
  element->get_vertex_indices(ev_i);
  short w, v, i, j, sum_w, mat_index, k;

  int fv_ind=0; // index in array free_vtces .

  // loop over all vertices in element.
  for (v=0; v<nve; ++v) {
    
    // finds out whether vertex v in the element is fixed or free,
    // as according to argument free_vtces[]
    bool free_vertex = false;
    for (k=0; k<num_free_vtx; ++k) {
      if ( ev_i[v] == pd.get_vertex_index(free_vtces[k]) )
        free_vertex = true;
    }

    // If vertex is fixed, enters null blocks for that column.
    // Note that null blocks for the row will be taken care of by
    // the gradient null entries. 
    if (free_vertex==false) {
      for (w=0; w<nve; ++w) {
        if (v>=w) {
          sum_w = w*(w+1)/2; // 1+2+3+...+w
          mat_index = w*nve+v-sum_w;
          hessian[mat_index] = 0.;
        }
      }
    }
    else  {
    // If vertex is free, use finite difference on the gradient to find the Hessian.
      for (j=0;j<3;++j) {
        // perturb the coordinates of the vertex v in the j direction by delta
        vj_coord = (*free_vtces[fv_ind])[j];
        (*free_vtces[fv_ind])[j]+=delta;
        //compute the gradient at the perturbed point location
        valid = this->compute_element_gradient_expanded(pd, element, free_vtces,
                              grad_vec1, num_free_vtx, metric_value, err); MSQ_CHKERR(err);
        assert(valid);
        //compute the numerical Hessian
        for (w=0; w<nve; ++w) {
          if (v>=w) {
            //finite difference to get some entries of the Hessian
            fd = (grad_vec1[w]-grad_vec[w])/delta;
            // For the block at position w,v in a matrix, we need the corresponding index
            // (mat_index) in a 1D array containing only upper triangular blocks.
            sum_w = w*(w+1)/2; // 1+2+3+...+w
            mat_index = w*nve+v-sum_w;
          
            for (i=0; i<3; ++i)
              hessian[mat_index][i][j] = fd[i];   
     
          }
        }
        // put the coordinates back where they belong
        (*free_vtces[fv_ind])[j] = vj_coord;
      }
      ++fv_ind;
    }
  }

  delete[] grad_vec1;

  FUNCTION_TIMER_END();
  return true;
}


#undef __FUNC__
#define __FUNC__ "QualityMetric::compute_vertex_numerical_gradient"
/*!  Numerically calculates the gradient of a vertex-based QualityMetric
  value on the given free vertex.  The metric is evaluated at MsqVertex
  'vertex', and the gradient is calculated with respect to the degrees
  of freedom associated with MsqVertices in the 'vertices' array.
*/
bool QualityMetric::compute_vertex_numerical_gradient(PatchData &pd,
                                                      MsqVertex &vertex,
                                                      MsqVertex* free_vtces[],
                                                      Vector3D grad_vec[],
                                                      int num_vtx,
                                                      double &metric_value,
                                                      MsqError &err)
{
   /*!TODO: (MICHAEL)  Try to inline this function (currenlty conflicts
      with MsqVertex.hpp).*/    
  MSQ_DEBUG_PRINT(2,"Computing Gradient (QualityMetric's numeric, vertex based.\n");
  
  bool valid=this->evaluate_vertex(pd, &(vertex), metric_value, err);
  MSQ_CHKERR(err);

  if (!valid)
    return false;
  
  double delta = 10e-6;
  double metric_value1=0;
  for (int v=0; v<num_vtx; ++v) 
  {
    /* gradient in the x, y, z direction */
    for (int j=0;j<3;++j) 
    {
      // perturb the coordinates of the free vertex in the j direction by delta
      (*free_vtces[v])[j]+=delta;
      //compute the function at the perturbed point location
      this->evaluate_vertex(pd, &(vertex),  metric_value1, err); MSQ_CHKERR(err);
      //compute the numerical gradient
      grad_vec[v][j]=(metric_value1-metric_value)/delta;
      // put the coordinates back where they belong
      (*free_vtces[v])[j] -= delta;
    }
  }
  return true;  
}

