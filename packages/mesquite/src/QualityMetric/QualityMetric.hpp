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

    (2006) kraftche@cae.wisc.edu    
   
  ***************************************************************** */

/*! \file QualityMetric.hpp
    \brief
Header file for the Mesquite::QualityMetric class

  \author Thomas Leurent
  \author Michael Brewer
  \date   2002-05-01
 */

#ifndef QualityMetric_hpp
#define QualityMetric_hpp

#ifndef MSQ_USE_OLD_C_HEADERS
#  include <cmath>
#else
#  include <math.h>
#endif

#ifndef MSQ_USE_OLD_STD_HEADERS
#  include <vector>
#  include <algorithm>
#else
#  include <vector.h>
#  include <algorithm.h>
#endif 

#include "Mesquite.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"

#ifdef _MSC_VER
   typedef unsigned uint32_t;
#elif defined(HAVE_STDINT_H)
#  include <stdint.h>
#elif defined(HAVE_INTTYPES_H)
#  include <inttypes.h>
#endif

namespace Mesquite
{
   
     /*! \class QualityMetric
       \brief Base class for concrete quality metrics.
     */
   class PatchData;
   class MsqMeshEntity;
   
   MESQUITE_EXPORT class QualityMetric
   {
   protected:

     QualityMetric( ) {}

   public:

     enum MetricType
     {
        VERTEX_BASED,  /**< Iterate over vertices to evaluate metric. */
        ELEMENT_BASED  /**< Iterate over elements to evaluate metric. */
     };

     virtual ~QualityMetric()
      {}
     
     virtual MetricType get_metric_type() const = 0;
     
     virtual msq_std::string get_name() const = 0;

      //! 1 if metric should be minimized, -1 if metric should be maximized.
     virtual int get_negate_flag() const = 0;
     
      /**\brief Get locations at which metric can be evaluated
       *
       * Different metrics are evaluated for different things within
       * a patch.  For example, an element-based metric will be evaluated
       * once for each element in patch, a vertex-based metric once for 
       * each veretx, and a target/sample-point based metric will be 
       * evaluated once for each samle point in each element.  This method
       * returns a list of handles, one for each location in the patch
       * at which the metric can be evaluated.  The handle values are used
       * as input to the evaluate methods.
       *\param pd       The patch
       *\param handles  Output list of handles
       *\param free_vertices_only If true, only pass back evaluation points
       *         that depend on at least one free vertex.
       */
     virtual
     void get_evaluations( PatchData& pd, 
                           msq_std::vector<size_t>& handles, 
                           bool free_vertices_only,
                           MsqError& err ) = 0;
     
     /**\brief Get metric value at a logical location in the patch.
      *
      * Evaluate the metric at one location in the PatchData.
      *\param pd     The patch.
      *\param handle The location in the patch (as passed back from get_evaluations).
      *\param value  The output metric value.
      */
     virtual
     bool evaluate( PatchData& pd, 
                    size_t handle, 
                    double& value, 
                    MsqError& err ) = 0;
     
     /**\brief Get metric value at a logical location in the patch.
      *
      * Evaluate the metric at one location in the PatchData.
      *\param pd      The patch.
      *\param handle  The location in the patch (as passed back from get_evaluations).
      *\param value   The output metric value.
      *\param indices The free vertices that the evaluation is a function
      *               of, specified as vertex indices in the PatchData.
      */
     virtual
     bool evaluate_with_indices( PatchData& pd,
                    size_t handle,
                    double& value,
                    msq_std::vector<size_t>& indices,
                    MsqError& err ) = 0;
     
     /**\brief Get metric value and gradient at a logical location in the patch.
      *
      * Evaluate the metric at one location in the PatchData.
      *\param pd      The patch.
      *\param handle  The location in the patch (as passed back from get_evaluations).
      *\param value   The output metric value.
      *\param indices The free vertices that the evaluation is a function
      *               of, specified as vertex indices in the PatchData.
      *\param gradient The gradient of the metric as a function of the
      *               coordinates of the free vertices passed back in
      *               the indices list.
      */
     virtual
     bool evaluate_with_gradient( PatchData& pd,
                    size_t handle,
                    double& value,
                    msq_std::vector<size_t>& indices,
                    msq_std::vector<Vector3D>& gradient,
                    MsqError& err );
     
     /**\brief Get metric value and gradient at a logical location in the patch.
      *
      * Evaluate the metric at one location in the PatchData.
      *\param pd      The patch.
      *\param handle  The location in the patch (as passed back from get_evaluations).
      *\param value   The output metric value.
      *\param indices The free vertices that the evaluation is a function
      *               of, specified as vertex indices in the PatchData.
      *\param gradient The gradient of the metric as a function of the
      *               coordinates of the free vertices passed back in
      *               the indices list.
      *\param Hessian_diagonal The 3x3 blocks along the diagonal of
      *               the Hessian matrix.
      */
     virtual
     bool evaluate_with_Hessian_diagonal( PatchData& pd,
                    size_t handle,
                    double& value,
                    msq_std::vector<size_t>& indices,
                    msq_std::vector<Vector3D>& gradient,
                    msq_std::vector<SymMatrix3D>& Hessian_diagonal,
                    MsqError& err );
     
     /**\brief Get metric value and deravitives at a logical location in the patch.
      *
      * Evaluate the metric at one location in the PatchData.
      *\param pd      The patch.
      *\param handle  The location in the patch (as passed back from get_evaluations).
      *\param value   The output metric value.
      *\param indices The free vertices that the evaluation is a function
      *               of, specified as vertex indices in the PatchData.
      *\param gradient The gradient of the metric as a function of the
      *               coordinates of the free vertices passed back in
      *               the indices list.
      *\param Hessian The Hessian of the metric as a function of the 
      *               coordinates. The Hessian is passed back as the
      *               upper-triangular portion of the matrix in row-major
      *               order, where each Matrix3D is the portion of the
      *               Hessian with respect to the vertices at the
      *               corresponding positions in the indices list.
      */
     virtual
     bool evaluate_with_Hessian( PatchData& pd,
                    size_t handle,
                    double& value,
                    msq_std::vector<size_t>& indices,
                    msq_std::vector<Vector3D>& gradient,
                    msq_std::vector<Matrix3D>& Hessian,
                    MsqError& err );

       //!Escobar Barrier Function for Shape and Other Metrics
       // det = signed determinant of Jacobian Matrix at a Vertex
       // delta = scaling parameter
     static inline double vertex_barrier_function(double det, double delta) 
            { return 0.5*(det+sqrt(det*det+4*delta*delta)); }
  //protected:

      /** \brief Remove from vector any gradient terms corresponding 
       *         to a fixed vertex.
       *
       * Remove terms from vector that correspond to fixed vertices.
       *\param type            Element type
       *\param fixed_vertices  Bit flags, one per vertex, 1 if
       *                       vertex is fixed.
       *\param gradients       Array of gradients
       */
      static void remove_fixed_gradients( EntityTopology type, 
                                          uint32_t fixed_vertices, 
                                          msq_std::vector<Vector3D>& gradients );

      /** \brief Remove from vectors any gradient terms and hessian
       *         diagonal blcoks corresponding to a fixed vertex.
       *
       * Remove terms from vector that correspond to fixed vertices.
       *\param type            Element type
       *\param fixed_vertices  Bit flags, one per vertex, 1 if
       *                       vertex is fixed.
       *\param gradients       Array of gradients
       *\param hess_diagonal_blocks   Array of diagonal blocks of Hessian matrix.
       */
      static void remove_fixed_diagonals( EntityTopology type, 
                                          uint32_t fixed_vertices, 
                                          msq_std::vector<Vector3D>& gradients,
                                          msq_std::vector<SymMatrix3D>& hess_diagonal_blocks );

      /** \brief Remove from vector any Hessian blocks corresponding 
       *         to a fixed vertex.
       *
       * Remove blocks from vector that correspond to fixed vertices.
       *\param type            Element type
       *\param fixed_vertices  Bit flags, one per vertex, 1 if
       *                       vertex is fixed.
       *\param hessians        Array of Hessian blocks (upper trianguler, row-major)
       */
      static void remove_fixed_hessians ( EntityTopology type, 
                                          uint32_t fixed_vertices, 
                                          msq_std::vector<Matrix3D>& hessians );
     
     /** \brief Convert fixed vertex format from list to bit flags
      *
      * Given list of pointers to fixed vertices as passed to
      * evaluation functions, convert to bit flag format used
      * for many utility functions in this class.  Bits correspond
      * to vertices in the canonical vertex order, beginning with
      * the least-significant bit.  The bit is cleared for free
      * vertices and set (1) for fixed vertices.
      */
      static uint32_t fixed_vertex_bitmap( PatchData& pd, 
                                           const MsqMeshEntity* elem,
                                           msq_std::vector<size_t>& free_indices );
      
     
     //! takes an array of coefficients and an array of metrics (both of length num_value)
     //! and averages the contents using averaging method 'method'.
      double weighted_average_metrics(const double coef[],
                                    const double metric_values[],
                                    const int& num_values, MsqError &err);

       /*!AveragingMethod allows you to set how the quality metric values
         attained at each sample point will be averaged together to produce
         a single metric value for an element.
       */
     enum AveragingMethod
     {
        LINEAR,                 //!< the linear average
        RMS,                    //!< the root-mean-squared average
        HMS,                    //!< the harmonic-mean-squared average
        SUM,                    //!< the sum of the values
        SUM_SQUARED,            //!< the sum of the squares of the values
        HARMONIC,               //!< the harmonic average
        LAST_WITH_HESSIAN=HARMONIC,
        MINIMUM,                //!< the minimum value
        MAXIMUM,                //!< the maximum value
        GEOMETRIC,              //!< the geometric average
        LAST_WITH_GRADIENT=GEOMETRIC,
        STANDARD_DEVIATION,     //!< the standard deviation squared of the values
        MAX_OVER_MIN,           //!< the maximum value minus the minum value
        MAX_MINUS_MIN,          //!< the maximum value divided by the minimum value
        SUM_OF_RATIOS_SQUARED   //!< (1/(N^2))*(SUM (SUM (v_i/v_j)^2))
     };

  private:
     int feasible;
     
     msq_std::vector<Matrix3D> tmpHess;
   };


} //namespace


#endif // QualityMetric_hpp
