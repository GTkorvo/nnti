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

/*! \file DistanceFromTarget.hpp

Header file for the Mesquite::DistanceFromTarget class

  \author Thomas Leurent
  \date   2004-09-29
 */


#ifndef DistanceFromTarget_hpp
#define DistanceFromTarget_hpp

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"
#include "ElementQM.hpp"
#include "MsqMatrix.hpp"
#include "TargetCalculator.hpp"
#include "WeightCalculator.hpp"

namespace Mesquite
{
  /*! \class DistanceFromTarget
    \brief Base class for the computation of the distance from target between
           the target matrices W and the actual corner matrices A. 
  */
  class DistanceFromTarget : public ElementQM
  {
  public:
  
    DistanceFromTarget( TargetCalculator* tc,
                        WeightCalculator* wc );
    
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~DistanceFromTarget();

    static const SamplePoints* get_dft_sample_pts();
    
    virtual msq_std::string get_name() const;
    void set_name( const msq_std::string& s ) { mName = s; }
    virtual int get_negate_flag() const;

  protected:
    void get_W_matrices( size_t elem_index, PatchData& pd,
                         Matrix3D T[], size_t num_T,
                         double c_k[], MsqError& err );
                                             
      //! For a given element, compute each corner matrix A, and given a target
      //! corner matrix W, returns \f$ T=AW^{-1} \f$ for each corner.    
    void compute_T_matrices( size_t elem_index, PatchData& pd,
                           Matrix3D T[], size_t num_T, double c_k[], MsqError &err);
 
    bool get_barrier_function(PatchData& pd, const double &tau, double &h, MsqError &err);

  private:
    
    TargetCalculator* targetCalc;
    WeightCalculator* weightCalc;
    msq_std::string mName;
  };
  
  inline DistanceFromTarget::DistanceFromTarget( TargetCalculator* tc,
                                                 WeightCalculator* wc )
    : targetCalc(tc), weightCalc(wc)
    {}

  inline void DistanceFromTarget::get_W_matrices( size_t elem_index,
                                                  PatchData& pd,
                                                  Matrix3D T[], 
                                                  size_t num_T,
                                                  double c_k[],
                                                  MsqError& err )
  {
    const unsigned num_corner = pd.element_by_index(elem_index).corner_count();
    if (num_corner > num_T) {
      MSQ_SETERR(err)("Insufficient array length", MsqError::INVALID_ARG);
      return;
    }
    const SamplePoints* pts = get_dft_sample_pts();
    
    for (unsigned i = 0; i < num_corner; ++i)
    {
      targetCalc->get_3D_target( pd, elem_index, pts, i, 
                                 *reinterpret_cast<MsqMatrix<3,3>*>(T+i), 
                                 err );
      MSQ_ERRRTN(err);
      
      c_k[i] = weightCalc->get_weight( pd, elem_index, pts, i, err );
      MSQ_ERRRTN(err);
    }
  }
    
  
  
  inline void DistanceFromTarget::compute_T_matrices(size_t elem_index, PatchData& pd,
                        Matrix3D T[], size_t num_T, double c_k[], MsqError &err)
  {    
    const unsigned num_corner = pd.element_by_index(elem_index).corner_count();
    pd.element_by_index(elem_index).compute_corner_matrices( pd, T, num_T, err );
    MSQ_ERRRTN(err);
    const SamplePoints* pts = get_dft_sample_pts();
    
    MsqMatrix<3,3> W;
    for (unsigned i = 0; i < num_corner; ++i)
    {
      targetCalc->get_3D_target( pd, elem_index, pts, i, W, err );
      MSQ_ERRRTN(err);
      timesInvA( T[i], *reinterpret_cast<Matrix3D*>(&W) );
      
      c_k[i] = weightCalc->get_weight( pd, elem_index, pts, i, err );
      MSQ_ERRRTN(err);
    }
  }


  /*! Returns the 
   */
  inline bool DistanceFromTarget::get_barrier_function(PatchData& pd, const double &tau, double &h, MsqError &err)
  { 

     double delta=pd.get_barrier_delta(err);  MSQ_ERRZERO(err);

     // Note: technically, we want delta=eta*tau-max
     //       whereas the function above gives delta=eta*alpha-max
     //      
     //       Because the only requirement on eta is eta << 1,
     //       and because tau-max = alpha-max/0.707 we can
     //       ignore the discrepancy

     if (delta==0) { 
        if (tau < MSQ_DBL_MIN ) {
           return false;
        }
        else {
           h=tau;
        }

     // Note: when delta=0, the vertex_barrier_function
     //       formally gives h=tau as well.
     //       We just do it this way to avoid any 
     //       roundoff issues.
     // Also: when delta=0, this metric is identical
     //       to the original condition number with
     //       the barrier at tau=0

     }
     else {
        h = 0.5*(tau+sqrt(tau*tau+4*delta*delta));

        if (h<MSQ_DBL_MIN && fabs(tau) > MSQ_DBL_MIN ) { 
          h = delta*delta/fabs(tau); }
 
        // Note: Analytically, h is strictly positive, but
        //       it can be zero numerically if tau
        //       is a large negative number 
        //       In the h=0 case, we use a different analytic
        //       approximation to compute h.
     }
     if (h<MSQ_DBL_MIN) {
       MSQ_SETERR(err)("Barrier function is zero due to excessively large "
                       "negative area compared to delta.\nTry to untangle "
                       "mesh another way.", MsqError::INVALID_MESH);
       return false;
     }
     return true;
  }  
  
} //namespace


#endif // DistanceFromTarget_hpp
