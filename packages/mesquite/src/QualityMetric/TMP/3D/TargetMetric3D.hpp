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


/** \file TargetMetric3D.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TARGET_METRIC_3D_HPP
#define MSQ_TARGET_METRIC_3D_HPP

#include "Mesquite.hpp"
#include <string>

namespace MESQUITE_NS {

class MsqError;
template <unsigned R, unsigned C> class MsqMatrix;

/**\brief A metric for comparing a 3x3 matrix A with a 3x3 target matrix W
 *
 * Implement a scalar function \f$\mu(A,W)\f$ where A and W are 3x3 matrices.
 */
class TargetMetric3D {
public:
    // used by code templatized to work with either this class or TargetMetric2D
  enum { MATRIX_DIM = 3 };

  MESQUITE_EXPORT virtual
  ~TargetMetric3D();

  MESQUITE_EXPORT virtual
  std::string get_name() const = 0;

    /**\brief Evaluate \f$\mu(A,W)\f$
     *
     *\param A 3x3 active matrix
     *\param W 3x3 target matrix
     *\param result Output: value of function
     *\return false if function cannot be evaluated for given A and W
     *          (e.g. division by zero, etc.), true otherwise.
     */
  MESQUITE_EXPORT virtual
  bool evaluate( const MsqMatrix<3,3>& A, 
                 const MsqMatrix<3,3>& W, 
                 double& result, 
                 MsqError& err ) = 0;
  
    /**\brief Gradient of \f$\mu(A,W)\f$ with respect to components of A
     *
     *\param A 3x3 active matrix
     *\param W 3x3 target matrix
     *\param result Output: value of function
     *\param deriv_wrt_A Output: partial deriviatve of \f$\mu\f$ wrt each term of A,
     *                           evaluated at passed A.
     *                           \f[\left[\begin{array}{ccc} 
     *                            \frac{\partial\mu}{\partial A_{0,0}} & 
     *                            \frac{\partial\mu}{\partial A_{0,1}} & 
     *                            \frac{\partial\mu}{\partial A_{0,2}} \\ 
     *                            \frac{\partial\mu}{\partial A_{1,0}} & 
     *                            \frac{\partial\mu}{\partial A_{1,1}} & 
     *                            \frac{\partial\mu}{\partial A_{1,2}} \\ 
     *                            \frac{\partial\mu}{\partial A_{2,0}} & 
     *                            \frac{\partial\mu}{\partial A_{2,1}} & 
     *                            \frac{\partial\mu}{\partial A_{2,2}}
     *                            \end{array}\right]\f]
     *\return false if function cannot be evaluated for given A and W
     *          (e.g. division by zero, etc.), true otherwise.
     */
  MESQUITE_EXPORT virtual
  bool evaluate_with_grad( const MsqMatrix<3,3>& A,
                           const MsqMatrix<3,3>& W,
                           double& result,
                           MsqMatrix<3,3>& deriv_wrt_A,
                           MsqError& err );
  
    /**\brief Hessian of \f$\mu(A,W)\f$ with respect to components of A
     *
     *\param A 3x3 active matrix
     *\param W 3x3 target matrix
     *\param result Output: value of function
     *\param deriv_wrt_A Output: partial deriviatve of \f$\mu\f$ wrt each term of A,
     *                           evaluated at passed A.
     *\param second_wrt_A Output: 9x9 matrix of second partial deriviatve of \f$\mu\f$ wrt 
     *                           each term of A, in row-major order.  The symmetric 
     *                           matrix is decomposed into 3x3 blocks and only the upper diagonal
     *                           blocks, in row-major order, are returned.
     *                           \f[\left[\begin{array}{ccc|ccc|ccc}
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}^2} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}\partial A_{0,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}\partial A_{0,2}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}\partial A_{1,0}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}\partial A_{1,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}\partial A_{1,2}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}\partial A_{2,0}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}\partial A_{2,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}\partial A_{2,2}} \\
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}\partial A_{0,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,1}^2} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,1}\partial A_{0,2}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,1}\partial A_{1,0}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,1}\partial A_{1,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,1}\partial A_{1,2}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,1}\partial A_{2,0}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,1}\partial A_{2,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,1}\partial A_{2,2}} \\
     *                           \frac{\partial^{2}\mu}{\partial A_{0,0}\partial A_{0,2}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,1}\partial A_{0,2}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,2}^2} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,2}\partial A_{1,0}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,2}\partial A_{1,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,2}\partial A_{1,2}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,2}\partial A_{2,0}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,2}\partial A_{2,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{0,2}\partial A_{2,2}} \\
     *                           \hline & & &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,0}^2} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,0}\partial A_{1,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,0}\partial A_{1,2}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,0}\partial A_{2,0}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,0}\partial A_{2,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,0}\partial A_{2,2}} \\
     *                           & & &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,0}\partial A_{1,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,1}^2} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,1}\partial A_{1,2}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,1}\partial A_{2,0}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,1}\partial A_{2,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,1}\partial A_{2,2}} \\
     *                           & & &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,0}\partial A_{1,2}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,1}\partial A_{1,2}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,2}^2} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,2}\partial A_{2,0}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,2}\partial A_{2,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{1,2}\partial A_{2,2}} \\
     *                           \hline & & & & & &
     *                           \frac{\partial^{2}\mu}{\partial A_{2,0}^2} &
     *                           \frac{\partial^{2}\mu}{\partial A_{2,0}\partial A_{2,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{2,0}\partial A_{2,2}} \\
     *                           & & & & & &
     *                           \frac{\partial^{2}\mu}{\partial A_{2,0}\partial A_{2,1}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{2,1}^2} &
     *                           \frac{\partial^{2}\mu}{\partial A_{2,1}\partial A_{2,2}} \\
     *                           & & & & & &
     *                           \frac{\partial^{2}\mu}{\partial A_{2,0}\partial A_{2,2}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{2,1}\partial A_{2,2}} &
     *                           \frac{\partial^{2}\mu}{\partial A_{2,2}^2} \\
     *                            \end{array}\right]\f]
     *        
     *\return false if function cannot be evaluated for given A and W
     *          (e.g. division by zero, etc.), true otherwise.
     */
  MESQUITE_EXPORT virtual
  bool evaluate_with_hess( const MsqMatrix<3,3>& A,
                           const MsqMatrix<3,3>& W,
                           double& result,
                           MsqMatrix<3,3>& deriv_wrt_A,
                           MsqMatrix<3,3> second_wrt_A[6],
                           MsqError& err );
                           
protected:
  static inline bool invalid_determinant( double d )
    { return d < 1e-12; }
};

} // namespace Mesquite

#endif
