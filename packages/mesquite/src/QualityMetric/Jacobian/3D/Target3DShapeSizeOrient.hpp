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


/** \file Target3DShapeSizeOrient.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TARGET_3D_SHAPE_SIZE_ORIENT_HPP
#define MSQ_TARGET_3D_SHAPE_SIZE_ORIENT_HPP

#include "Mesquite.hpp"
#include "TargetMetric3D.hpp"

namespace Mesquite {

/**\brief |T-I|^2
 *
 * A target metric for volume elements that optimizes
 * element shape, size, and orientation.
 */
class Target3DShapeSizeOrient : public TargetMetric3D
{
public:

  MESQUITE_EXPORT virtual
  ~Target3DShapeSizeOrient();

  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<3,3>& A, 
                 const MsqMatrix<3,3>& W, 
                 double& result, 
                 MsqError& err );
  
};


} // namespace Mesquite

#endif
