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


/** \file DistanceFromTarget.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "DistanceFromTarget.hpp"
#include "SamplePoints.hpp"

namespace Mesquite {

const SamplePoints* DistanceFromTarget::get_dft_sample_pts()
{ 
  static const SamplePoints samplePoints( true, false, false, false );
  return &samplePoints; 
}

msq_std::string DistanceFromTarget::get_name() const { return mName; }
int DistanceFromTarget::get_negate_flag() const { return 1; }
DistanceFromTarget::~DistanceFromTarget() {}

} // namespace Mesquite
