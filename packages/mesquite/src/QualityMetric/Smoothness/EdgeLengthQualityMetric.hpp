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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file EdgeLengthQualityMetric.hpp

Header file for the Mesquite::EdgeLengthQualityMetric class

  \author Michael Brewer
  \date   2002-06-13
 */


#ifndef EdgeLengthQualityMetric_hpp
#define EdgeLengthQualityMetric_hpp
#include "MsqMeshEntity.hpp"
#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "SmoothnessQualityMetric.hpp"
#include "Vector3D.hpp"

namespace Mesquite
{
     /*! \class EdgeLengthQualityMetric
       \brief Computes the lengths of the edges connected to given a vertex..
       
        EdgeLengthQualityMetric is a vertex based metric which computes
        the lengths of the edges connected to a given vertex and then
        averages those together, using the specified averaging method
        The metric uses SUM as the default averaging method.
     */
   class MsqMeshEntity;
   class MsqVertex;
   
   class EdgeLengthQualityMetric : public SmoothnessQualityMetric
  {
   public:
    
    EdgeLengthQualityMetric()
       {
         avgMethod=SUM;
         feasible=0;
         set_metric_type(QualityMetric::VERTEX_BASED);
         set_name("Edge Length Metric");
       }
    
      // virtual destructor ensures use of polymorphism during destruction
    virtual ~EdgeLengthQualityMetric()
       {}

  
   protected:
   
    bool evaluate_vertex(PatchData &pd, MsqVertex *vert, double &fval,
                         MsqError &err);
  };


} //namespace


#endif // EdgeLengthQualityMetric_hpp


