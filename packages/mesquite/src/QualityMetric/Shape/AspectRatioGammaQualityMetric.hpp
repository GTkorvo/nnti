// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file AspectRatioGammaQualityMetric.hpp
  \brief
  Header file for the Mesquite::AspectRatioGammaQualityMetric class

  \author Michael Brewer
  \date   2002-05-16
 */


#ifndef AspectRatioGammaQualityMetric_hpp
#define AspectRatioGammaQualityMetric_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ShapeQualityMetric.hpp"
#include "Vector3D.hpp"
#ifdef USE_C_PREFIX_INCLUDES
#include <cmath>
#else
#include <math.h>
#endif
namespace Mesquite
{
     /*! \class AspectRatioGammaQualityMetric
       \brief Object for computing the aspect ratio gamma of
       simplicial elements.
     */
   class AspectRatioGammaQualityMetric : public ShapeQualityMetric
   {
   public:
       /*! The function create_new is used to create a shape quality metric.
        */
     static ShapeQualityMetric* create_new()
        {
          ShapeQualityMetric* m = new AspectRatioGammaQualityMetric();
          return m;
        }
     
       //! virtual destructor ensures use of polymorphism during destruction
     virtual ~AspectRatioGammaQualityMetric()
        {}
     
   protected:
     
     
   private:
       //constants used in metric calculations
     double fourDivRootThree;
     double twelveDivRootTwo;
     AspectRatioGammaQualityMetric()
        {
          fourDivRootThree=4.0/sqrt(3.0);
          twelveDivRootTwo=12.0/sqrt(2.0);
          feasible=0;
          evalMode=QualityMetric::ELEMENT_VERTICES;
          set_metric_type(QualityMetric::ELEMENT_BASED);
          set_name("Aspect Ratio Gamma");
        }
       //!Computes the aspect ratio gamma of element.  If element
       //!is not a tetrahedron or triangle, sets an error.
     bool evaluate_element(PatchData& pd,
                           MsqMeshEntity* element, double &fval,
                           MsqError &err);
   };
   
   
} //namespace


#endif // AspectRatioGammaQualityMetric_hpp


