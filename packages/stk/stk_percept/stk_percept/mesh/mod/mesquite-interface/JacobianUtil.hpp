/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef JacobianUtil_hpp
#define JacobianUtil_hpp

#if !defined(__IBMCPP__) && defined(STK_BUILT_IN_SIERRA)

#include <Mesquite.hpp>
#include <MsqError.hpp>
#include <MsqDebug.hpp>
#include <AveragingQM.hpp>
#include <Vector3D.hpp>
#include <MsqMatrix.hpp>
#include <Exponent.hpp>

#include <stk_percept/PerceptMesh.hpp>

namespace stk {
  namespace percept {

    using namespace Mesquite;

    class JacobianUtil : public AveragingQM
    {
      const double a2Con;
      const Exponent b2Con;
      const Exponent c2Con;
      
      const double a3Con;
      const Exponent b3Con;
      const Exponent c3Con;
      Vector3D mCoords[4];

    public:
      double   mMetrics[8];
      MsqMatrix<3,3> mJ[8];

      JacobianUtil() : 
        AveragingQM( QualityMetric::LINEAR ),

        a2Con(1.0), 
        b2Con(0.0), 
        c2Con(1.0),
        a3Con(1.0),
        b3Con(0.0),
        c3Con(1.0)
      {
      }

      bool operator()(double& averageJ, PerceptMesh& eMesh, stk::mesh::Entity& element, stk::mesh::FieldBase *coord_field=0);

    };
  }
}

#endif
#endif
