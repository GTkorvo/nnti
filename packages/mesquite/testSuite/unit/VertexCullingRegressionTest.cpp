// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Michael Brewer
//       ORG: Sandia National Labs
//    E-MAIL: mbrewer@sandia.gov
//
// ORIG-DATE: May 8, 2003
//  LAST-MOD: 23-Jul-03 at 17:44:57 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file VertexCullingRegressionTest.cpp

Regression testing using the vertex culling algorithms. 
 */
// DESCRIP-END.
//

#include "PatchDataInstances.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"
#include <math.h>

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "QualityAssessor.hpp"
#include "MsqMessage.hpp"

#include "EdgeLengthQualityMetric.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "LaplacianSmoother.hpp"
#include "PlanarDomain.hpp"
#include "TerminationCriterion.hpp"
#include "MeshImpl.hpp"
using namespace Mesquite;

class VertexCullingRegressionTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(VertexCullingRegressionTest);
  CPPUNIT_TEST (test_laplacian_smoothing_with_cull);
  CPPUNIT_TEST_SUITE_END();
  
private:
  int pF;//PRINT_FLAG
public:
  void setUp()
  {
      //pF=1;//PRINT_FLAG IS ON
      pF=0;//PRINT_FLAG IS OFF
  }

  void tearDown()
  {
  }
  
public:
  VertexCullingRegressionTest()
    {}

  void test_laplacian_smoothing_with_cull()
    {
        /* Read a VTK Mesh file */
      MsqError err;
      Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
      mesh->read_vtk("../../meshFiles/2D/VTK/square_quad_10_rand.vtk", err);
        //Make sure no errors
      CPPUNIT_ASSERT(!err.errorOn);
        // initialises a MeshSet object
      MeshSet mesh_set1;
      mesh_set1.add_mesh(mesh, err); MSQ_CHKERR(err);
      
        // creates an intruction queue
      InstructionQueue queue1;
      
        // creates a mean ratio quality metric ...
      ShapeQualityMetric* shape_metric = new ConditionNumberQualityMetric;
      SmoothnessQualityMetric* lapl_met = new EdgeLengthQualityMetric;
      lapl_met->set_averaging_method(QualityMetric::RMS,err);
      
        // creates the laplacian smoother  procedures
      LaplacianSmoother* lapl1 = new LaplacianSmoother(err);
      LaplacianSmoother* lapl2 = new LaplacianSmoother(err);
      QualityAssessor stop_qa=QualityAssessor(shape_metric,QualityAssessor::MAXIMUM);
      stop_qa.add_quality_assessment(lapl_met,QualityAssessor::ALL_MEASURES,err);
      
        //**************Set termination criterion****************
      TerminationCriterion sc2;
      sc2.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1000,err);
        //Make sure no errors
      CPPUNIT_ASSERT(!err.errorOn);
      sc2.add_criterion_type_with_double(TerminationCriterion::SUCCESSIVE_IMPROVEMENTS_ABSOLUTE,0.0,err);
        //Make sure no errors
      CPPUNIT_ASSERT(!err.errorOn);
        //set a criterion with a culling method for the inner criterion
      TerminationCriterion sc_cull;
      sc_cull.set_culling_type(TerminationCriterion::VERTEX_MOVEMENT_ABSOLUTE,.1,err);
      TerminationCriterion sc_cull_2;
      sc_cull_2.set_culling_type(TerminationCriterion::VERTEX_MOVEMENT_ABSOLUTE,.000001,err);
        //Make sure no errors
      CPPUNIT_ASSERT(!err.errorOn);
      lapl1->set_outer_termination_criterion(&sc2);
      lapl2->set_outer_termination_criterion(&sc2);
      lapl1->set_inner_termination_criterion(&sc_cull);
      lapl2->set_inner_termination_criterion(&sc_cull_2);
        // sets a culling method on the first QualityImprover
      lapl1->add_culling_method(PatchData::NO_BOUNDARY_VTX);
      lapl2->add_culling_method(PatchData::NO_BOUNDARY_VTX);
        // adds 1 pass of pass1 to mesh_set1
      queue1.add_quality_assessor(&stop_qa,err);
       //Make sure no errors
      CPPUNIT_ASSERT(!err.errorOn);
      queue1.add_preconditioner(lapl1,err);
        //Make sure no errors
      CPPUNIT_ASSERT(!err.errorOn);
      queue1.set_master_quality_improver(lapl2, err);
       //Make sure no errors
      CPPUNIT_ASSERT(!err.errorOn);
      queue1.add_quality_assessor(&stop_qa,err);
        //Make sure no errors
      CPPUNIT_ASSERT(!err.errorOn);
      queue1.run_instructions(mesh_set1, err);
      CPPUNIT_ASSERT(!err.errorOn);
      
        //Make sure no errors
      CPPUNIT_ASSERT(!err.errorOn);
      delete shape_metric;
      delete lapl1;
      delete lapl2;
      delete lapl_met;
    }
  
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(VertexCullingRegressionTest, "VertexCullingRegressionTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(VertexCullingRegressionTest, "Regression");
