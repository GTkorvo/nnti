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
  \file   MultiplyQualityMetric.cpp
  \brief 

  \author Michael Brewer
  \date   2002-05-09
*/

#include "MultiplyQualityMetric.hpp"
#include "QualityMetric.hpp"
#include "Vector3D.hpp"
#include "MsqMessage.hpp"

using namespace Mesquite;
#undef __FUNC__
#define __FUNC__ "MultiplyQualityMetric::MultiplyQualityMetric"

MultiplyQualityMetric::MultiplyQualityMetric(QualityMetric* qm1, QualityMetric* qm2, MsqError &err){
  if(qm1 == NULL || qm2 == NULL){
    err.set_msg("MultiplyQualityMetric constructor passed NULL pointer.");
  }
  set_qmetric2(qm2);
  set_qmetric1(qm1);
  feasible=qm1->get_feasible_constraint();
  if(qm2->get_feasible_constraint())
    feasible=qm2->get_feasible_constraint();
  int n_flag=qm1->get_negate_flag();
  if(n_flag!=qm2->get_negate_flag()){
    PRINT_WARNING("MultiplyQualityMetric is being used to compose a metric that should be minimized\n with a metric that should be maximized.");
    set_negate_flag(1);
  }
  else{
    set_negate_flag(n_flag);
  }

  // Checks that metrics are of the same type 
  if ( qm1->get_metric_type() != qm2->get_metric_type() ) {
    err.set_msg("Cannot multiply a vertex-based QM with an element-based QM.");
  } else {
    set_metric_type(qm1->get_metric_type());
  }
  
  set_name("Composite Multiply");
}

#undef __FUNC__
#define __FUNC__ "MultiplyQualityMetric::evaluate_element"

/*! Returns qMetric1->evaluate_element(element, err) multiplied by
  qMetric2-evaluate_element(element, err)*/
bool MultiplyQualityMetric::evaluate_element(PatchData& pd,
                                             MsqMeshEntity *element,
                                             double &value,
                                             MsqError &err)
{
  bool valid_flag;
  double metric1, metric2;  
  valid_flag=get_qmetric1()->evaluate_element(pd, element, metric1, err);
  if(!valid_flag)
    return false;
  valid_flag=get_qmetric2()->evaluate_element(pd, element, metric2, err); 
  value = metric1*metric2;
    //if the first metric was invalid we have already returned
    //so we return whatever the flag was on the second metric.
  return
    valid_flag;
}

#undef __FUNC__
#define __FUNC__ "MultiplyQualityMetric::evaluate_vertex"
/*! Returns qMetric1->evaluate_vertex(...) multiplied by
  qMetric2-evaluate_vertex(...)*/
bool MultiplyQualityMetric::evaluate_vertex(PatchData& pd,
                                            MsqVertex* vert,
                                            double &value,
                                            MsqError& err)
{
  bool valid_flag;
  double metric1, metric2;  
  valid_flag=get_qmetric1()->evaluate_vertex(pd, vert, metric1, err);
  if(!valid_flag)
    return false;
  valid_flag=get_qmetric2()->evaluate_vertex(pd, vert, metric2, err);
  value = metric1*metric2;
    //if the first metric was invalid we have already returned
    //so we return whatever the flag was on the second metric.
  return
    valid_flag;
}



