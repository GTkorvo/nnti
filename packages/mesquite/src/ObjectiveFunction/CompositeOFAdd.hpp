// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file CompositeOFAdd.hpp

Header file for the Mesquite:: CompositeOFAdd class

  \author Michael Brewer
  \date   2002-06-24
 */


#ifndef CompositeOFAdd_hpp
#define CompositeOFAdd_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "ObjectiveFunction.hpp"
#include "PatchData.hpp"
#include <list>

namespace Mesquite
{
     /*!\class CompositeOFAdd
       \brief Adds two ObjectiveFunction values together.
     */
   class MsqMeshEntity;
   class PatchData;
   class CompositeOFAdd : public ObjectiveFunction
   {
   public:
     CompositeOFAdd(ObjectiveFunction*, ObjectiveFunction*);
     virtual ~CompositeOFAdd();
     virtual bool concrete_evaluate(PatchData &patch, double &fval,
                                    MsqError &err);
     virtual std::list<QualityMetric*> get_quality_metric_list();
     
	protected:
     bool compute_analytical_gradient(PatchData &patch,Vector3D *const &grad,
                                      double &OF_val,  MsqError &err,
				      size_t array_size);
	private:
     ObjectiveFunction* objFunc1;
     ObjectiveFunction* objFunc2;
   };
}//namespace
#endif //  CompositeOFAdd_hpp
