/*!
  \file   MeanRatioQualityMetric.cpp
  \brief  

  \author Michael Brewer
  \date   2002-06-9
*/
#include <vector>
#include "MeanRatioQualityMetric.hpp"
#include <math.h>
#include "Vector3D.hpp"
#include "ShapeQualityMetric.hpp"
#include "QualityMetric.hpp"

using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "MeanRatioQualityMetric::MeanRatioQualityMetric"

MeanRatioQualityMetric::MeanRatioQualityMetric()
{
  set_negate_flag(1);
  avgMethod=QualityMetric::LINEAR;
  feasible=1;
  evalMode=QualityMetric::ELEMENT_VERTICES;
  set_name("Mean Ratio");
}

bool MeanRatioQualityMetric::evaluate_element(PatchData &pd,
                                              MsqMeshEntity *element,
                                              double &fval,
                                              MsqError &err)
{
  double metric_values[20];
  fval=0.0;
  bool return_flag;
  std::vector<size_t> v_i;
  element->get_vertex_indices(v_i);
    //only 3 temp_vec will be sent to mean ratio calculator, but the
    //additional vector3D may be needed during the calculations
  Vector3D temp_vec[5];
  MsqVertex *vertices=pd.get_vertex_array(err);
  switch(element->get_element_type()){
    case TRIANGLE:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[2]=vertices[v_i[2]]-vertices[v_i[0]];
        //make relative to equilateral
      temp_vec[1]=((2*temp_vec[2])-temp_vec[0])*MSQ_SQRT_THREE_INV;
      return_flag=mean_ratio_2d(temp_vec,fval,err);
      return return_flag;
      break;
    case QUADRILATERAL:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[1]=vertices[v_i[3]]-vertices[v_i[0]];
      return_flag=mean_ratio_2d(temp_vec,metric_values[0],err);
      if(!return_flag)
        return false;
      
      temp_vec[0]=vertices[v_i[2]]-vertices[v_i[1]];
      temp_vec[1]=vertices[v_i[0]]-vertices[v_i[1]];
      return_flag=mean_ratio_2d(temp_vec,metric_values[1],err);
      if(!return_flag)
        return false;
      temp_vec[0]=vertices[v_i[3]]-vertices[v_i[2]];
      temp_vec[1]=vertices[v_i[1]]-vertices[v_i[2]];
      return_flag=mean_ratio_2d(temp_vec,metric_values[2],err);
      if(!return_flag)
        return false;
      temp_vec[0]=vertices[v_i[0]]-vertices[v_i[3]];
      temp_vec[1]=vertices[v_i[2]]-vertices[v_i[3]];
      return_flag=mean_ratio_2d(temp_vec,metric_values[3],err);
      if(!return_flag)
        return false;
      fval=average_metrics(metric_values,4,err);
      break;
    case TETRAHEDRON:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[3]=vertices[v_i[2]]-vertices[v_i[0]];
      temp_vec[4]=vertices[v_i[3]]-vertices[v_i[0]];
        //transform to equilateral tet
      temp_vec[1]=((2*temp_vec[3])-temp_vec[0])/MSQ_SQRT_THREE;
      temp_vec[2]=((3*temp_vec[4])-temp_vec[0]-temp_vec[3])/
        (MSQ_SQRT_THREE*MSQ_SQRT_TWO);
      return_flag=mean_ratio_3d(temp_vec,fval,err);
      break;
    case HEXAHEDRON:
      temp_vec[0]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[1]=vertices[v_i[3]]-vertices[v_i[0]];
      temp_vec[2]=vertices[v_i[4]]-vertices[v_i[0]];
      return_flag=mean_ratio_3d(temp_vec,metric_values[0],err);
      if(!return_flag)
        return false;
      temp_vec[0]=vertices[v_i[2]]-vertices[v_i[1]];
      temp_vec[1]=vertices[v_i[0]]-vertices[v_i[1]];
      temp_vec[2]=vertices[v_i[5]]-vertices[v_i[1]];
      return_flag=mean_ratio_3d(temp_vec,metric_values[1],err);
      if(!return_flag)
        return false;
      temp_vec[0]=vertices[v_i[3]]-vertices[v_i[2]];
      temp_vec[1]=vertices[v_i[1]]-vertices[v_i[2]];
      temp_vec[2]=vertices[v_i[6]]-vertices[v_i[2]];
      return_flag=mean_ratio_3d(temp_vec,metric_values[2],err);
      if(!return_flag)
        return false;
      temp_vec[0]=vertices[v_i[0]]-vertices[v_i[3]];
      temp_vec[1]=vertices[v_i[2]]-vertices[v_i[3]];
      temp_vec[2]=vertices[v_i[7]]-vertices[v_i[3]];
      return_flag=mean_ratio_3d(temp_vec,metric_values[3],err);
      if(!return_flag)
        return false;
      temp_vec[0]=vertices[v_i[7]]-vertices[v_i[4]];
      temp_vec[1]=vertices[v_i[5]]-vertices[v_i[4]];
      temp_vec[2]=vertices[v_i[0]]-vertices[v_i[4]];
      return_flag=mean_ratio_3d(temp_vec,metric_values[4],err);
      if(!return_flag)
        return false;
      temp_vec[0]=vertices[v_i[4]]-vertices[v_i[5]];
      temp_vec[1]=vertices[v_i[6]]-vertices[v_i[5]];
      temp_vec[2]=vertices[v_i[1]]-vertices[v_i[5]];
      return_flag=mean_ratio_3d(temp_vec,metric_values[5],err);
      if(!return_flag)
        return false;
      temp_vec[0]=vertices[v_i[5]]-vertices[v_i[6]];
      temp_vec[1]=vertices[v_i[7]]-vertices[v_i[6]];
      temp_vec[2]=vertices[v_i[2]]-vertices[v_i[6]];
      return_flag=mean_ratio_3d(temp_vec,metric_values[6],err);
      if(!return_flag)
        return false;
      temp_vec[0]=vertices[v_i[6]]-vertices[v_i[7]];
      temp_vec[1]=vertices[v_i[4]]-vertices[v_i[7]];
      temp_vec[2]=vertices[v_i[3]]-vertices[v_i[7]];
      return_flag=mean_ratio_3d(temp_vec,metric_values[0],err);
      if(!return_flag)
        return false;
      fval=average_metrics(metric_values,8,err);
      
      break;
    default:
      fval=0.0;
  }// end switch over element type
  return true;
}


