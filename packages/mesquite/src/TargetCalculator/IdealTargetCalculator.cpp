/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file IdealTargetCalculator.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "IdealTargetCalculator.hpp"
#include "MappingFunction.hpp"
#include "PatchData.hpp"
#include "MsqError.hpp"
#include "IdealElements.hpp"
#include "ElemSampleQM.hpp"

namespace Mesquite {

bool IdealTargetCalculator::surface_targets_are_3D() const
  { return false; }
 
bool IdealTargetCalculator:: get_3D_target( PatchData& pd, 
                                            size_t element,
                                            const SamplePoints* pts,
                                            unsigned sample,
                                            MsqMatrix<3,3>& W,
                                            MsqError& err )
{
  MsqMeshEntity& elem = pd.element_by_index( element );
  EntityTopology type = elem.get_element_type();
  unsigned dim = ElemSampleQM::side_dim_from_sample( sample );
  unsigned num = ElemSampleQM::side_num_from_sample( sample );
  const MappingFunction3D* func = pd.get_mapping_function_3D( type );

  const Vector3D* verts = unit_element( type );
  if (!verts) {
      MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT);
      return false;
  }
  
  jc.get_Jacobian_3D( func, 0, dim, num, verts, elem.node_count(), W, err );
  MSQ_ERRZERO(err);
  return true;
}

bool IdealTargetCalculator:: get_2D_target( PatchData& pd, 
                                            size_t element,
                                            const SamplePoints* pts,
                                            unsigned sample,
                                            MsqMatrix<3,2>& W,
                                            MsqError& err )
{
  MsqMeshEntity& elem = pd.element_by_index( element );
  EntityTopology type = elem.get_element_type();
  unsigned dim = ElemSampleQM::side_dim_from_sample( sample);
  unsigned num = ElemSampleQM::side_num_from_sample( sample );
  const MappingFunction2D* func = pd.get_mapping_function_2D( type );

  const Vector3D* verts = unit_element( type );
  if (!verts) {
      MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT);
      return false;
  }
  
  jc.get_Jacobian_2D( func, 0, dim, num, verts, elem.node_count(), W, err );
  MSQ_ERRZERO(err);
  
  if (orientSurfElems) {
    Vector3D n;
    switch (dim) {
      case 0: pd.get_domain_normal_at_corner( element, num, n, err ); break;
      case 1: pd.get_domain_normal_at_mid_edge( element, num, n, err ); break;
      default:pd.get_domain_normal_at_element( element, n, err ); break;
    }
    
      // Assume ideal element is in XY plane
    MsqMatrix<3,3> rotation;
    n.normalize();
    const double pl = n[0]*n[0] + n[1]*n[1];
    if (pl < DBL_EPSILON) {
      if (n[2] < 0.0) {
        W(0,0) = -W(0,0);
        W(0,1) = -W(0,1);
      }
      return true;
    }
   
    const double f = (1.0 - n[2])/pl;
    rotation(0,0) = n[2] + f * n[1] * n[1];
    rotation(1,1) = n[2] + f * n[0] * n[0];
    rotation(1,0) = rotation(0,1) = -f * n[0] * n[1];
    rotation(0,2) =  n[0];
    rotation(1,2) =  n[1];
    rotation(2,2) =  n[2];
    rotation(2,1) = -n[1];
    rotation(2,0) = -n[0];
    
    W = rotation * W;
  }
  
  return true;
}

} // namespace Mesquite
