/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file Target2DUntangleMu.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TARGET_2D_UNTANGLE_MU_HPP
#define MSQ_TARGET_2D_UNTANGLE_MU_HPP

#include "Mesquite.hpp"
#include "TargetMetric2D.hpp"

namespace MESQUITE_NS {

/**\brief Composite untangle metric
 *
 * This metric should be combined with Target2DSize or Target2DShapeSize
 * to produce a concrete untangle metric.
 *
 * \f$ \mu^\prime(T) = (|d| - d)^2 \f$
 * \f$ d(T) = \sigma - \epsilon - \mu(T() \f$
 *
 */
class Target2DUntangleMu : public TargetMetric2D
{
private:
  TargetMetric2D* mBaseMetric;
  double mConstant;

public:

  Target2DUntangleMu( TargetMetric2D* base, 
                      double sigma = 1.0 ) 
    : mBaseMetric(base),
      mConstant(0.99*sigma) /* default epsilon is 0.01*sigma */
    {}

  Target2DUntangleMu( TargetMetric2D* base, 
                      double sigma,
                      double epsilon ) 
    : mBaseMetric(base),
      mConstant(sigma-epsilon) 
    {}

  MESQUITE_EXPORT virtual
  ~Target2DUntangleMu();

  MESQUITE_EXPORT virtual
  std::string get_name() const;

  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<2,2>& A, 
                 const MsqMatrix<2,2>& W, 
                 double& result, 
                 MsqError& err );
   
 MESQUITE_EXPORT virtual
  bool evaluate_with_grad( const MsqMatrix<2,2>& A,
                           const MsqMatrix<2,2>& W,
                           double& result,
                           MsqMatrix<2,2>& deriv_wrt_A,
                           MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_hess( const MsqMatrix<2,2>& A,
                           const MsqMatrix<2,2>& W,
                           double& result,
                           MsqMatrix<2,2>& deriv_wrt_A,
                           MsqMatrix<2,2> second_wrt_A[3],
                           MsqError& err );
};


} // namespace MESQUITE_NS

#endif
