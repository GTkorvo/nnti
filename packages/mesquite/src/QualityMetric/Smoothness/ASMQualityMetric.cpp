/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
/*!
  \file   ASMQualityMetric.cpp
  \brief  QualityMetric class for ASM (Area Smoothness Metric)

  \author Michael Brewer
  \date   2002-06-9
*/
#include <vector>
#include "ASMQualityMetric.hpp"
#include <math.h>
#include "Vector3D.hpp"
#include "ShapeQualityMetric.hpp"
#include "QualityMetric.hpp"
#include "MsqMessage.hpp"
using namespace Mesquite;
MSQ_USE(vector);

#undef __FUNC__
#define __FUNC__ "ASMQualityMetric::ASMQualityMetric"

ASMQualityMetric::ASMQualityMetric()
{
  MsqError err;
  set_metric_type(ELEMENT_BASED);
  set_element_evaluation_mode(ELEMENT_VERTICES, err);
  MSQ_CHKERR(err);
  avgMethod=QualityMetric::MAXIMUM;
  feasible=0;
  set_name("Area Smoothness");
}

bool ASMQualityMetric::evaluate_element(PatchData &pd,
                                                    MsqMeshEntity *element,
                                                    double &fval,
                                                    MsqError &err)
{
  double temp_double;
  size_t elem_ind=pd.get_element_index(element);
  vector<size_t> adj_elems;
 
  MsqMeshEntity *elems = pd.get_element_array(err);
  switch(element->get_element_type()){
    case TRIANGLE:
    case QUADRILATERAL:
      pd.get_adjacent_entities_via_n_dim(1,elem_ind,adj_elems,err);
      break;
    case TETRAHEDRON:
    case HEXAHEDRON:
      pd.get_adjacent_entities_via_n_dim(2,elem_ind,adj_elems,err);
      break;
    default:
      err.set_msg("ASM quality metric not implemented for this element type.");
  };
  int num_samp=adj_elems.size();
  if(num_samp < 1){
    fval=0.0;
  }
  else{
    double* met_vals = new double [num_samp];
    int i=0;
    switch(element->get_element_type()){
      case TRIANGLE:
      case QUADRILATERAL:
        temp_double=element->compute_unsigned_area(pd,err);
          //PRINT_INFO("\nunsigned area = %f",temp_double);
        for(i=0;i<num_samp;++i){
          met_vals[i]=elems[adj_elems[i]].compute_unsigned_area(pd,err);
            //PRINT_INFO("neighboring nunsigned area = %f",met_vals[i]);
          if((temp_double+met_vals[i])>MSQ_MIN){
            met_vals[i]=fabs((temp_double-met_vals[i])/(temp_double+
                                                        met_vals[i]));
          }
          else
            met_vals[i]=0.0;
        }
        break;                                
        
      case TETRAHEDRON:
      case HEXAHEDRON:
        temp_double=element->compute_unsigned_volume(pd,err);
        for(i=0;i<num_samp;++i){
          met_vals[i]=elems[adj_elems[i]].compute_unsigned_volume(pd,err);
          if((temp_double+met_vals[i])>MSQ_MIN){
            met_vals[i]=fabs((temp_double-met_vals[i])/(temp_double+
                                                        met_vals[i]));
          }
          else
            met_vals[i]=0.0;
        }
        break;
      default:
        err.set_msg("ASM quality metric not implemented for this element type.");
    };
    fval=average_metrics(met_vals,num_samp,err);
      //PRINT_INFO("\nRETURNING %f \n",fval);
    delete []met_vals;
  }
  return true;
}


