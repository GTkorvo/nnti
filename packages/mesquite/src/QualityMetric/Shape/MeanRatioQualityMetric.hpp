// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file MeanRatioQualityMetric.hpp

Header file for the Mesquite::MeanRatioQualityMetric class

\author Michael Brewer
\author Thomas Leurent
\date   2002-06-19
 */


#ifndef MeanRatioQualityMetric_hpp
#define MeanRatioQualityMetric_hpp
#include "MsqMeshEntity.hpp"
#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ShapeQualityMetric.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"
#include "PatchData.hpp"
//Michael delete
#include "MsqMessage.hpp"

namespace Mesquite
{
   /*! \class MeanRatioQualityMetric
     \brief Computes the mean ratio of given element.

     The metric does not use the sample point functionality or the
     compute_weighted_jacobian.  It evaluates the metric at
     the element vertices, and uses the isotropic ideal element.
     Optionally, the metric computation can be raised to the
     'pow_dbl' power.  This does not necessarily raise the metric
     value to the 'pow_dbl' power but instead raises each local
     metric.  For example, if the corner mean ratios of a quadraliteral
     element were m1,m2,m3, and m4 and we set pow_dbl=2 and
     used linear averaging, the metric value would then be
     m = .25(m1*m1 + m2*m2 + m3*m3 + m4*m4).  The metric does
     require a feasible region, and the metric needs to be minimized
     if pow_dbl is greater than zero and maximized if pow_dbl
     is less than zero.  pow_dbl being equal to zero is invalid.
   */
   class MeanRatioQualityMetric : public ShapeQualityMetric
   {
   public:
      MeanRatioQualityMetric(double pow_dbl=1.0) : ShapeQualityMetric() {
       MsqError err;
       set_metric_type(ELEMENT_BASED);
       set_element_evaluation_mode(ELEMENT_VERTICES, err); MSQ_CHKERR(err);

       set_gradient_type(ANALYTICAL_GRADIENT);
       set_hessian_type(ANALYTICAL_HESSIAN);
       avgMethod=QualityMetric::LINEAR;
       feasible=1;
       set_name("Mean Ratio");
       
         //Note:  the following are redundant since set_metric_power is called
       set_negate_flag(1);

       a2Con =  1.0 / 2.0;
       b2Con =  1.0;
       c2Con = -1.0;

       a3Con =  1.0 / 3.0;
       b3Con =  1.0;
       c3Con = -2.0 / 3.0;
         //the above are redundant since set_metric_power is called
       
       set_metric_power(pow_dbl);
      }

      //! virtual destructor ensures use of polymorphism during destruction
      virtual ~MeanRatioQualityMetric() {
      }
     
      //! evaluate using mesquite objects 
      bool evaluate_element(PatchData &pd, MsqMeshEntity *element, 
                            double &fval, MsqError &err); 

      bool compute_element_analytical_gradient(PatchData &pd,
                                               MsqMeshEntity *element,
                                               MsqVertex *free_vtces[], 
                                               Vector3D grad_vec[],
                                               int num_vtx, 
                                               double &metric_value,
                                               MsqError &err);

      bool compute_element_analytical_hessian(PatchData &pd,
                                              MsqMeshEntity *e,
                                              MsqVertex *v[], 
                                              Vector3D g[],
                                              Matrix3D h[],
                                              int nv, 
                                              double &m,
                                              MsqError &err);
      
    private:
       //! Sets the power value in the metric computation.
     void set_metric_power(double pow_dbl)
        {
          if(fabs(pow_dbl)<MSQ_MIN){
            pow_dbl=1.0;
            PRINT_WARNING("\nInvalid power passed to set_metric_power(double ), Metric using the default value, 1.0, instead.");
          }
          if(pow_dbl<0)
            set_negate_flag(-1);
          else
            set_negate_flag(1);
          a2Con=pow(.5,pow_dbl);
          b2Con=pow_dbl;
          c2Con=-pow_dbl;
          a3Con=pow(1.0/3.0,pow_dbl);
          b3Con=pow_dbl;
          c3Con=-2.0*pow_dbl/3.0;
        }
      // arrays used in Hessian computations 
      // We allocate them here, so that one allocation only is done.
      // This gives a big computation speed increase.
      Vector3D mCoords[4]; // Vertex coordinates for the (decomposed) elements
      Vector3D mGradients[32]; // Gradient of metric with respect to the coords
      Vector3D mAccumGrad[8];  // Accumulated gradients (composed merit)
      Matrix3D mHessians[80]; // Hessian of metric with respect to the coords
      double   mMetrics[8]; // Metric values for the (decomposed) elements
      double   g_factor[8]; // Metric values for the (decomposed) elements
     double   h_factor[8]; // Metric values for the (decomposed) elements
       //variables used in the definition of the metric (2d and 3d)
     double a2Con;
     double b2Con;
     double c2Con;
     
     double a3Con;
     double b3Con;
     double c3Con;
     
     
   };
} //namespace


#endif // MeanRatioQualityMetric_hpp
