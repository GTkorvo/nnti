// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 30-Oct-02 at 17:58:24 by Thomas Leurent
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#ifdef USE_STD_INCLUDES
#include <iostream>
#else
#include <iostream.h>
#endif

#ifdef USE_C_PREFIX_INCLUDES
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#include "MsqMessage.hpp"
#include "MsqTimer.hpp"
#include "Mesquite.hpp"
#include "TSTT_Base.h"
#include "MesquiteUtilities.hpp" //  for writeShowMeMesh()
#include "MesquiteError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorythms
#include "Randomize.hpp"
#include "MeanRatioQualityMetric.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "LPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"
#include "ConjugateGradient.hpp"
using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "main"
int main()
{     
  char file_name[128];
  // Reads a TSTT Mesh file 
  TSTT::Mesh_Handle mesh;
  TSTT::MeshError tstt_err;
  TSTT::Mesh_Create(&mesh, &tstt_err);
  strcpy(file_name, "../../meshFiles/2D/VTK/tangled_quad.vtk");
  
  TSTT::Mesh_Load(mesh, file_name, &tstt_err);
  
    // Mesquite error object
  MsqError err;
  
    // initialises a MeshSet object
  MeshSet mesh_set1;
  mesh_set1.add_mesh(mesh, err); MSQ_CHKERR(err);
  
    // creates an intruction queue
  InstructionQueue queue1;
  
    // creates a mean ratio quality metric ...
  ShapeQualityMetric* shape_metric = ConditionNumberQualityMetric::create_new();
  UntangleQualityMetric* untangle = UntangleBetaQualityMetric::create_new(2);
  Randomize* pass0 = new Randomize(.05);
    // ... and builds an objective function with it
    //LInfTemplate* obj_func = new LInfTemplate(shape_metric);
  LInfTemplate* obj_func = new LInfTemplate(untangle);
  LPTemplate* obj_func2 = new LPTemplate(shape_metric, 2, err);
    // creates the steepest descent optimization procedures
  obj_func2->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
  ConjugateGradient* pass1 = new ConjugateGradient( obj_func, err );
  
    //SteepestDescent* pass2 = new SteepestDescent( obj_func2 );
  ConjugateGradient* pass2 = new ConjugateGradient( obj_func2, err );
  pass2->set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH,err,1,1);
  pass2->set_patch_type(PatchData::GLOBAL_PATCH,err);
  QualityAssessor stop_qa=QualityAssessor(shape_metric,QualityAssessor::MAXIMUM);
  QualityAssessor stop_qa2=QualityAssessor(shape_metric,QualityAssessor::MAXIMUM);
  stop_qa2.add_quality_assessment(shape_metric,QualityAssessor::AVERAGE,err);
  
  stop_qa.add_quality_assessment(untangle,QualityAssessor::ALL_MEASURES,err);
  stop_qa.set_stopping_assessment(untangle,QualityAssessor::MAXIMUM,err);
    // **************Set stopping criterion**************
    //untangle beta should be 0 when untangled
  TerminationCriterion sc1;
  sc1.add_criterion_type_with_double(TerminationCriterion::QUALITY_IMPROVEMENT,
                                     .000001,err);
  TerminationCriterion sc3;
  sc3.add_criterion_type_with_int(TerminationCriterion::ITERATION_BOUND,10,
                                  err);
  TerminationCriterion sc_rand;
  sc_rand.add_criterion_type_with_int(TerminationCriterion::ITERATION_BOUND,1,
                                      err);

    //StoppingCriterion sc1(&stop_qa,-1.0,.0000001);
    //StoppingCriterion sc3(&stop_qa2,.9,1.00000001);
    //StoppingCriterion sc2(StoppingCriterion::NUMBER_OF_PASSES,10);
    //StoppingCriterion sc_rand(StoppingCriterion::NUMBER_OF_PASSES,1);
    //either until untangled or 10 iterations
  pass0->set_outer_termination_criterion(&sc_rand);
  pass1->set_outer_termination_criterion(&sc1);
  pass2->set_inner_termination_criterion(&sc3);
  
    // sets a culling method on the first QualityImprover
  pass0->add_culling_method(PatchData::NO_BOUNDARY_VTX);
  pass1->add_culling_method(PatchData::NO_BOUNDARY_VTX);
  pass2->add_culling_method(PatchData::NO_BOUNDARY_VTX);
    // adds 1 pass of pass1 to mesh_set1
  queue1.add_quality_assessor(&stop_qa,err); MSQ_CHKERR(err);
    //queue1.add_preconditioner(pass0,err);MSQ_CHKERR(err);
    //queue1.add_preconditioner(pass1,err);MSQ_CHKERR(err);
    //queue1.set_master_quality_improver(pass2, err); MSQ_CHKERR(err);
  queue1.set_master_quality_improver(pass1, err); MSQ_CHKERR(err);
  queue1.add_quality_assessor(&stop_qa2,err); MSQ_CHKERR(err);
  writeVtkMesh("original_mesh", mesh, err); MSQ_CHKERR(err);
  
    // launches optimization on mesh_set1
  queue1.run_instructions(mesh_set1, err); MSQ_CHKERR(err);
  
  writeVtkMesh("smoothed_mesh", mesh, err); MSQ_CHKERR(err);
  PRINT_TIMING_DIAGNOSTICS();
}
