/*!
  \file   MinTemplate.cpp
  \brief  

  This Objective Function is the minimum of the quality metrics
  total=min (x)
  \author Lori Freitag
  \date   2002-07-18
*/
#include <math.h>
#include "MinTemplate.hpp"
using  namespace Mesquite;  

#undef __FUNC__
#define __FUNC__ "MinTemplate::MinTemplate"

MinTemplate::MinTemplate(QualityMetric *qualitymetric){
   set_quality_metric(qualitymetric);
   set_feasible(qualitymetric->get_feasible_constraint());
}

#undef __FUNC__
#define __FUNC__ "MinTemplate::~MinTemplate"

//Lori:  need to clean up here
MinTemplate::~MinTemplate(){

}

#undef __FUNC__
#define __FUNC__ "MinTemplate::concrete_evaluate"

double MinTemplate::concrete_evaluate(PatchData &patch, MsqError &err){

  //Total value of objective function
  double total_value=0;
  double temp_value=0;

  //For elements in Patch
  int index;
  QualityMetric* currentQM = get_quality_metric();
  if(currentQM->get_metric_type()==QualityMetric::ELEMENT_BASED){
    
    int num_elements=patch.num_elements();
    MsqMeshEntity* elems=patch.get_element_array(err);
 
    for (index=0; index<num_elements; index++){

      //evaluate metric for this elem
      currentQM->evaluate_element(patch, &elems[index], temp_value, err);
      MSQ_CHKERR(err);

      if(temp_value<total_value)
        total_value=temp_value;

    }//end loop over elements
  }//end if not VERTEX

  else if (currentQM->get_metric_type()==QualityMetric::VERTEX_BASED) {

    int num_vertices=patch.num_vertices();
    MsqVertex* vertices=patch.get_vertex_array(err);
 
    for (index=0; index<num_vertices;index++){

      //evaluate metric for this vertex
      currentQM->evaluate_vertex(patch, &vertices[index], temp_value, err);
      MSQ_CHKERR(err);

      if(temp_value<total_value)
        total_value=temp_value;

    }//end loop over vertices
  }//end else VERTEX
  else {
    err.set_msg("Make sure MetricType is initialised in concrete QualityMetric constructor.");
  }
  
  return total_value;
}
	
	
