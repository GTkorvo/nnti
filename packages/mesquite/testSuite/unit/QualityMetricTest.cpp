// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Michael Brewer
//       ORG: Sandia National Labs
//    E-MAIL: mbrewer@sandia.gov
//
// ORIG-DATE: 03-Dec-02
//  LAST-MOD: 23-Jul-03 at 17:34:39 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file QualityMetricTest.cpp

Unit testing of various QualityMetrics primarily to test for
correct metric return values. 

\author Michael Brewer
\author Thomas Leurent

 */
// DESCRIP-END.
//

#include "Mesquite.hpp"
#include "PatchData.hpp"
#include "PatchDataInstances.hpp"
#include "QualityMetric.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "VertexConditionNumberQualityMetric.hpp"
#include "GeneralizedConditionNumberQualityMetric.hpp"
#include "MeanRatioQualityMetric.hpp"
#include "InverseMeanRatioQualityMetric.hpp"
#include "AspectRatioGammaQualityMetric.hpp"
#include "MultiplyQualityMetric.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "MsqMessage.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"
#include "ASMQualityMetric.hpp"
#include "SmoothnessQualityMetric.hpp"
#include "LocalSizeQualityMetric.hpp"
#include "CornerJacobianQualityMetric.hpp"
#include "PowerQualityMetric.hpp"
#include "ScalarAddQualityMetric.hpp"
#include "MeshImpl.hpp"
#include "MeshSet.hpp"
#include "PlanarDomain.hpp"

#include <math.h>

using namespace Mesquite;

using std::cout;
using std::cerr;
using std::endl;

class QualityMetricTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(QualityMetricTest);
    //Test condition number and generalized condition number metrics
  CPPUNIT_TEST (test_condition_number);
    //Test mean ratio and inverse mean ratio metrics
  CPPUNIT_TEST (test_mean_ratio);
    //Test apect ratio gamma (Tri's and Tet's only)
  CPPUNIT_TEST (test_aspect_ratio_gamma);
    //Test composite multiply
  CPPUNIT_TEST (test_composite_multiply);
  CPPUNIT_TEST (test_other_composites);
  
    //Test edge length metric
  CPPUNIT_TEST (test_edge_length_metric);  
    //Test averaging methods
  CPPUNIT_TEST (test_averaging_method);
    // test analytical gradients
  CPPUNIT_TEST (test_mean_ratio_tri_gradient_planar);
  CPPUNIT_WORK_IN_PROGRESS (test_mean_ratio_tri_gradient_nonplanar);
  CPPUNIT_TEST (test_mean_ratio_quad_gradient_planar);
  CPPUNIT_WORK_IN_PROGRESS (test_mean_ratio_quad_gradient_nonplanar);
  CPPUNIT_TEST (test_mean_ratio_tet_gradient);
  CPPUNIT_TEST (test_mean_ratio_hex_gradient);
    // test analytical Hessians
  CPPUNIT_TEST (test_mean_ratio_tri_hessian);
  CPPUNIT_TEST (test_mean_ratio_quad_hessian);
  CPPUNIT_TEST (test_mean_ratio_quad_hessian_linear);
  CPPUNIT_TEST (test_mean_ratio_quad_hessian_sum_squared);
  CPPUNIT_TEST (test_mean_ratio_quad_hessian_rms);
  CPPUNIT_TEST (test_mean_ratio_quad_hessian_harmonic);
  CPPUNIT_TEST (test_mean_ratio_quad_hessian_hms);
  CPPUNIT_TEST (test_mean_ratio_tet_hessian);
  CPPUNIT_TEST (test_mean_ratio_hex_hessian);
  CPPUNIT_TEST (test_mean_ratio_hex_hessian_linear);
  CPPUNIT_TEST (test_mean_ratio_hex_hessian_sum_squared);
  CPPUNIT_TEST (test_mean_ratio_hex_hessian_rms);
  CPPUNIT_TEST (test_mean_ratio_hex_hessian_harmonic);
  CPPUNIT_TEST (test_mean_ratio_hex_hessian_hms);
    //Test ASM (area smoothness quality metric)
  CPPUNIT_TEST (test_asm);
    //Test corner jacobian metric
  CPPUNIT_TEST (test_corner_jacobian_metric);
    //Test local size metric
  CPPUNIT_TEST (test_local_size_metric);
    //Test local size metric
  CPPUNIT_TEST (test_vertex_based_condition_number_metric);
  // test element gradient coming directly from the gradient
  // versus the one coming from the hessian
  CPPUNIT_TEST (test_mean_ratio_tri_grad_from_hessian);
  CPPUNIT_TEST (test_mean_ratio_quad_grad_from_hessian);
  CPPUNIT_TEST (test_mean_ratio_tet_grad_from_hessian);
  CPPUNIT_TEST (test_mean_ratio_hex_grad_from_hessian);

  CPPUNIT_TEST (test_untangle_metric);

  CPPUNIT_TEST_SUITE_END();
  
private:
  
  PatchData triPatch;
  PatchData quadPatch;
  PatchData tetPatch;
  PatchData hexPatch;
  PatchData invertedTri;
  PatchData invertedTet;
  PatchData idealTri;
  PatchData idealTet;
    //Tol used for double comparisons
  double qualTol;
  int pF;//PRINT_FLAG
public:
  void setUp()
  {
      //pF=1;//PRINT_FLAG IS ON
    pF=0;//PRINT_FLAG IS OFF
    MsqError err;
    
    qualTol = 1.e-12;
    
     /* Our triangular patch is made of two tris.  tri_1 is a perfect
        equilateral (the ideal for most metrics).  tri_2 is an arbitrary
        triangle.
     */
    create_qm_two_tri_patch_with_domain(triPatch, err);MSQ_CHKERR(err);
    
     /* Our quad patch is made of two quads.  quad_1 is a perfect
        square (the ideal for most metrics).  quad_2 is an arbitrary
        quad.
     */
    create_qm_two_quad_patch_with_domain(quadPatch,err);MSQ_CHKERR(err);
    
     /* Our tet patch is made of two tets.  tet_1 is a perfect
        equilateral (the ideal for most metrics).  tet_2 is an arbitrary
        tet.
     */
    create_qm_two_tet_patch(tetPatch,err);MSQ_CHKERR(err);
    
     /* Our hex patch is made of two hexes.  hex_1 is a perfect
        unit cube (the ideal for most metrics).  hex_2 is an arbitrary
        hex.
     */
     create_qm_two_hex_patch(hexPatch,err);MSQ_CHKERR(err);

       //'ideal' inverted tet
     create_one_inverted_tet_patch(invertedTet, err);MSQ_CHKERR(err);
       //ideal tri
     create_one_tri_patch(idealTri, err);MSQ_CHKERR(err);
       //ideal tet
     create_one_tet_patch(idealTet, err);MSQ_CHKERR(err);
  }

  void tearDown()
  {
    destroy_patch_with_domain(triPatch);
    destroy_patch_with_domain(quadPatch);
  }
  
public:
  QualityMetricTest()
    {}

#undef __FUNC__
#define __FUNC__ "QualityMetricTest::test_condition_number"
   void test_condition_number()
   {
     bool v_flag;
       //START WITH TRI's
     MsqError err;
     double val, val2;
     MsqMeshEntity* elems;
     MsqVertex* verts = triPatch.get_vertex_array(err);
     elems=triPatch.get_element_array(err);
     MSQ_CHKERR(err);
     ShapeQualityMetric *met = new ConditionNumberQualityMetric;
     ShapeQualityMetric *gmet = new GeneralizedConditionNumberQualityMetric;
       //Check condition number of ideal tri
     v_flag=met->evaluate_element(triPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
       //Check generalized condition number of ideal tri
     v_flag=gmet->evaluate_element(triPatch,&elems[0],val2,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val2,1.0,qualTol);
       //For now, make sure cond num and generalized cond num give
       //equivalent answer for arbitrary tri.
     v_flag=met->evaluate_element(triPatch,&elems[1],val,err); MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     v_flag=gmet->evaluate_element(triPatch,&elems[1],val2,err); MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     val -= val2;
     if(pF)
       PRINT_INFO("\nGEN TRI %f", val2);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,0.0,qualTol);
     
       //SECOND: QUAD's
     verts = quadPatch.get_vertex_array(err);
     elems = quadPatch.get_element_array(err);
       //Check condition number of ideal quad
     v_flag=met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
     
     //Check generalized condition number of ideal quad
     v_flag=gmet->evaluate_element(quadPatch,&elems[0],val2,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val2,1.0,qualTol);
       //For now, make sure cond num and generalized cond num give
       //equivalent answer for arbitrary quad.
     v_flag=met->evaluate_element(quadPatch,&elems[1],val,err); MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     v_flag=gmet->evaluate_element(quadPatch,&elems[1],val2,err); MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     
     val -= val2;
     if(pF)
       PRINT_INFO("\nGEN QUA %f", val2);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,0.0,qualTol);
     
       //THIRD TET's
     verts = tetPatch.get_vertex_array(err);
     elems = tetPatch.get_element_array(err);
       //Check condition number of ideal tet
     v_flag=met->evaluate_element(tetPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
       //Check generalized condition number of ideal tet
     v_flag=gmet->evaluate_element(tetPatch,&elems[0],val2,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val2,1.0,qualTol);
       //For now, make sure cond num and generalized cond num give
       //equivalent answer for arbitrary tet.
     v_flag=met->evaluate_element(tetPatch,&elems[1],val,err); MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     v_flag=gmet->evaluate_element(tetPatch,&elems[1],val2,err); MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     
     val -= val2;
     if(pF)
       PRINT_INFO("\nGEN TET %f", val2);
     
       //CPPUNIT_ASSERT(fabs(val)<qualTol);

       //FOURTH HEX's
     verts = hexPatch.get_vertex_array(err);
     elems = hexPatch.get_element_array(err);
       //Check condition number of ideal hex
     v_flag=met->evaluate_element(hexPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
       //Check generalized condition number of ideal hex
     v_flag=gmet->evaluate_element(hexPatch,&elems[0],val2,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val2,1.0,qualTol);
     
       //For now, make sure cond num and generalized cond num give
       //equivalent answer for arbitrary tet.
     v_flag=met->evaluate_element(hexPatch,&elems[1],val,err); MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     if(pF)
       PRINT_INFO("\nCON HEX %f", val);
     CPPUNIT_ASSERT(v_flag==true);
     
     v_flag=gmet->evaluate_element(hexPatch,&elems[1],val2,err); MSQ_CHKERR(err);
     CPPUNIT_ASSERT(v_flag==true);
     val -= val2;
     if(pF)
       PRINT_INFO("\nGEN HEX %f", val2);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,0.0,qualTol);
   }


   void test_mean_ratio()
   {
       //START WITH TRI's
     MsqError err;
     bool valid;
     double val;
     MsqMeshEntity* elems;
     MsqVertex* verts = triPatch.get_vertex_array(err);
     elems=triPatch.get_element_array(err);
     MSQ_CHKERR(err);
     ShapeQualityMetric *met = new MeanRatioQualityMetric;
     ShapeQualityMetric *imet = new InverseMeanRatioQualityMetric;
       //Check mean ratio of ideal tri
     met->evaluate_element(triPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nMEAN TRI %f", val);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
       //Check inverse mean ratio of ideal tri (INVERSE)
     imet->evaluate_element(triPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nInv MEAN TRI %f", val);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
       //SECOND: QUAD's
     verts = quadPatch.get_vertex_array(err);
     elems = quadPatch.get_element_array(err);
       //Check mean ratio of ideal quad
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nMEAN QUAD %f", val);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
       //Check inverse mean ratio of ideal quad (INVERSE)
     imet->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nInv MEAN QUAD %f", val);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);

       //THIRD TET's
     verts = tetPatch.get_vertex_array(err);
     elems = tetPatch.get_element_array(err);
       //Check mean ratio of ideal tet
     met->evaluate_element(tetPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nMEAN TET %f", val);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);

       //Check inverse mean ratio of ideal tet (INVERSE)
     imet->evaluate_element(tetPatch,&elems[0],val,err); MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nInv MEAN TET %f", val);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);

       //FOURTH HEX's
     verts = hexPatch.get_vertex_array(err);
     elems = hexPatch.get_element_array(err);
       //Check mean ratio of ideal hex
     valid = met->evaluate_element(hexPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(valid==true);
     if(pF)
       PRINT_INFO("\nMEAN HEX %f", val);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
       //Check inverse mean ratio of ideal hex (INVERSE)
     valid = imet->evaluate_element(hexPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(valid==true);
     if(pF)
       PRINT_INFO("\nInv MEAN HEX %f", val);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
   }
  
  void test_aspect_ratio_gamma()
   {
       //START WITH TRI's
     MsqError err;
     double val;
     bool valid=false;
     MsqMeshEntity* elems;
     MsqVertex* verts = triPatch.get_vertex_array(err);
     elems=triPatch.get_element_array(err);
     MSQ_CHKERR(err);
     ShapeQualityMetric *met = new AspectRatioGammaQualityMetric;
       //Check aspect ratio gamma of ideal tri
     valid=met->evaluate_element(triPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
     CPPUNIT_ASSERT(valid==true);
       //THIRD TET's
     verts = tetPatch.get_vertex_array(err);
     elems = tetPatch.get_element_array(err);
       //Check aspect ratio gamma of ideal tet
     valid=met->evaluate_element(tetPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT(valid==true);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
   }

  void test_composite_multiply()
   {
       //START WITH TRI's
     MsqError err;
     double val;
     MsqMeshEntity* elems;
     MsqVertex* verts = triPatch.get_vertex_array(err);
     elems=triPatch.get_element_array(err);
     MSQ_CHKERR(err);
     ShapeQualityMetric *mmet = new MeanRatioQualityMetric;
     ShapeQualityMetric *cmet = new ConditionNumberQualityMetric;
     CompositeQualityMetric *met = new MultiplyQualityMetric(mmet,
                                                            cmet,
                                                            err);
       //Check ideal tri
     met->evaluate_element(triPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nMULT TRI %f", val);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
       //SECOND: QUAD's
     verts = quadPatch.get_vertex_array(err);
     elems = quadPatch.get_element_array(err);
       //Check ideal quad
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nMULT QUAD %f", val);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
       //THIRD TET's
     verts = tetPatch.get_vertex_array(err);
     elems = tetPatch.get_element_array(err);
       //Check ideal tet
     
     if(!met->evaluate_element(tetPatch,&elems[0],val,err)){
        PRINT_INFO("\nMULT RETURNING FALSE FOR TET");
     }
     if(pF){
       PRINT_INFO("\nMULT TET %f", val);
       double tetm=0;
       if(!mmet->evaluate_element(tetPatch,&elems[0],tetm,err))
          PRINT_INFO("\nMMET TET %f", tetm);
       if(!cmet->evaluate_element(tetPatch,&elems[0],tetm,err))
          PRINT_INFO("\nCMET TET %f", tetm);
     }
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
       //FOURTH HEX's
     verts = hexPatch.get_vertex_array(err);
     elems = hexPatch.get_element_array(err);
       //Check ideal hex
     met->evaluate_element(hexPatch,&elems[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nMULT HEX %f", val);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
   }
  void test_other_composites()
     {
       MsqError err;
       double val;
       double temp_val;
       MsqMeshEntity* elems;
       MsqVertex* verts = triPatch.get_vertex_array(err);
       elems=triPatch.get_element_array(err);
       MSQ_CHKERR(err);
         //vertex based
       SmoothnessQualityMetric *m1 = new EdgeLengthQualityMetric;
       CompositeQualityMetric *pow_m1 = new PowerQualityMetric(m1,2,err);
       CompositeQualityMetric *sa_m1 = new ScalarAddQualityMetric(m1,2,err);
       m1->evaluate_vertex(triPatch,&verts[2],val,err);MSQ_CHKERR(err);
       pow_m1->evaluate_vertex(triPatch,&verts[2],temp_val,err);
       MSQ_CHKERR(err);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(val*val,temp_val,qualTol);
       sa_m1->evaluate_vertex(triPatch,&verts[2],temp_val,err);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(val+2.0,temp_val,qualTol);
         //element based
       ShapeQualityMetric *m2 = new ConditionNumberQualityMetric;
       CompositeQualityMetric *pow_m2 = new PowerQualityMetric(m2,2,err);
       CompositeQualityMetric *sa_m2 = new ScalarAddQualityMetric(m2,2,err);
       m2->evaluate_element(triPatch,&elems[0],val,err);MSQ_CHKERR(err);
       pow_m2->evaluate_element(triPatch,&elems[0],temp_val,err);
       MSQ_CHKERR(err);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(val*val,temp_val,qualTol);
       sa_m2->evaluate_element(triPatch,&elems[0],temp_val,err);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(val+2.0,temp_val,qualTol);
         //element based with a negative power
       CompositeQualityMetric *pow_mneg2 = new PowerQualityMetric(m2,-2,err);
       m2->evaluate_element(triPatch,&elems[0],val,err);MSQ_CHKERR(err);
       pow_mneg2->evaluate_element(triPatch,&elems[0],temp_val,err);
       MSQ_CHKERR(err);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/(val*val),temp_val,qualTol);
       delete m1;
       delete m2;
       delete pow_m1;
       delete pow_m2;
       delete sa_m1;
       delete sa_m2;
       delete pow_mneg2;
     }
  

  void test_edge_length_metric()
   {
       //START WITH TRI's
     MsqError err;
     double val;
     MsqMeshEntity* elems;
     MsqVertex* verts = triPatch.get_vertex_array(err);
     elems=triPatch.get_element_array(err);
     MSQ_CHKERR(err);
     SmoothnessQualityMetric *met = new EdgeLengthQualityMetric;
       //SmoothnessQualityMetric *lam_met = new EdgeLengthQualityMetric(EdgeLengthQualityMetric::DIVIDE_BY_LAMBDA);
       //Check aspect ratio gamma of ideal tri patch
       //Vert[2] has two edges connected of length 1
     met->evaluate_vertex(triPatch,&verts[2],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nEdge Length Metric tris (should be 2) %f", val);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,2.0,qualTol);
       //lam_met->evaluate_vertex(triPatch,&verts[2],val,err);MSQ_CHKERR(err);
       //if(pF)
       //PRINT_INFO("\nEdge Length Metric tris (should be 4) %f", val);
       //CPPUNIT_ASSERT_DOUBLES_EQUAL(val,4.0,qualTol);
       //Vert[1] has two edges connected of length 1
       //THIRD TET's
     verts = tetPatch.get_vertex_array(err);
     elems = tetPatch.get_element_array(err);
       //Check aspect ratio gamma of ideal tet
     met->evaluate_vertex(tetPatch,&verts[0],val,err);MSQ_CHKERR(err);
     if(pF)
       PRINT_INFO("\nEdge Length Metric tets (should be 3) %f", val);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,3.0,qualTol);
     delete met;
       //delete lam_met;
   }
  
  void test_asm()
     {
       QualityMetric *met = new ASMQualityMetric;
       MsqError err;
         //Test the metric on a single elemnt patch
       PatchData p1;
       create_one_tri_patch(p1, err);
       MsqMeshEntity* elem1=p1.get_element_array(err);
       double first_val;
       bool first_bool=met->evaluate_element(p1,&elem1[0],first_val,err);
         //Any non-connected element should have asm of 0.0
       CPPUNIT_ASSERT_DOUBLES_EQUAL(first_val, 0.0, qualTol);
       CPPUNIT_ASSERT(first_bool==true);
         //Test on a patch with two ideal tris
       PatchData p2;
       create_two_tri_patch(p2, err);
       MsqMeshEntity* elem2=p2.get_element_array(err);
       double second_val;
       bool second_bool=met->evaluate_element(p2,&elem2[0],second_val,err);
         //Two neighboring tris with equal area should have asm of 0.0
       CPPUNIT_ASSERT_DOUBLES_EQUAL(second_val, 0.0, qualTol);
       CPPUNIT_ASSERT(second_bool==true);

        //Test on a patch with two tris one not ideal
       PatchData p3;
       create_qm_two_tri_patch_with_domain(p3, err);
       MsqMeshEntity* elem3=p3.get_element_array(err);
       double third_val;
       bool third_bool=met->evaluate_element(p3,&elem3[0],third_val,err);
         //Two neighboring tris with equal area should have asm of 0.0
       CPPUNIT_ASSERT(third_val>0.0);
       CPPUNIT_ASSERT(third_val<1.0);
       CPPUNIT_ASSERT(third_bool==true);
       destroy_patch_with_domain(p3);
     }
  
  void test_corner_jacobian_metric()
     {
       QualityMetric *met = new CornerJacobianQualityMetric;
       MsqError err;
         //Test the metric on a single elemnt patch
       PatchData p1;
       create_one_tri_patch(p1, err);
       MsqMeshEntity* elem1=p1.get_element_array(err);
       double first_val;
       bool first_bool=met->evaluate_element(p1,&elem1[0],first_val,err);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(first_val, sqrt(3.0)/4.0, qualTol);
       CPPUNIT_ASSERT(first_bool==true);
         //Test on a patch with two ideal tris
       PatchData p2;
       create_two_tri_patch(p2, err);
       MsqMeshEntity* elem2=p2.get_element_array(err);
       double second_val;
       bool second_bool=met->evaluate_element(p2,&elem2[0],second_val,err);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(second_val, sqrt(3.0)/4.0, qualTol);
       CPPUNIT_ASSERT(second_bool==true);

        //Test on a patch with two tris one not ideal
       PatchData p3;
       create_qm_two_tri_patch_with_domain(p3, err);
       MsqMeshEntity* elem3=p3.get_element_array(err);
       double third_val;
       bool third_bool=met->evaluate_element(p3,&elem3[0],third_val,err);
       CPPUNIT_ASSERT(third_val>0.0);
       CPPUNIT_ASSERT(third_val<1.0);
       CPPUNIT_ASSERT(third_bool==true);
       destroy_patch_with_domain(p3);
     }
  
  void test_local_size_metric()
     {
       VolumeQualityMetric *met = new LocalSizeQualityMetric;
       MsqError err;
         //Test the metric on a single elemnt patch
       PatchData p1;
       create_one_tri_patch(p1, err);
       MsqVertex* vert1=p1.get_vertex_array(err);
       double first_val=-1;
       bool first_bool=met->evaluate_vertex(p1,&vert1[0],first_val,err);
       if(pF)
         PRINT_INFO("\nLocal size ideal tri %f", first_val);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(first_val, 1.0, qualTol);
       CPPUNIT_ASSERT(first_bool==true);
       first_bool=met->evaluate_vertex(p1,&vert1[1],first_val,err);
       if(pF)
         PRINT_INFO("\nLocal size ideal tri %f", first_val);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(first_val, 1.0, qualTol);
       CPPUNIT_ASSERT(first_bool==true);
         //Test on a patch with two ideal tris
       PatchData p2;
       create_two_tri_patch(p2, err);
       MsqVertex* vert2=p2.get_vertex_array(err);
       double second_val;
       bool second_bool=met->evaluate_vertex(p2,&vert2[0],second_val,err);
         //Two neighboring tris with equal area should have local size of 1.0
       if(pF)
         PRINT_INFO("\nLocal Size Metric two ideal tris %f", second_val);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(second_val, 1.0, qualTol);
       CPPUNIT_ASSERT(second_bool==true);
         /*
        //Test on a patch with two tris one not ideal
       PatchData p3;
       create_qm_two_tri_patch(p3, err);
       MsqMeshEntity* elem3=p3.get_element_array(err);
       double third_val;
       bool third_bool=met->evaluate_element(p3,&elem3[0],third_val,err);
         //Two neighboring tris with equal area should have asm of 0.0
       CPPUNIT_ASSERT(third_val>0.0);
       CPPUNIT_ASSERT(third_val<1.0);
       CPPUNIT_ASSERT(third_bool==true);
         */
     }
  void test_vertex_based_condition_number_metric()
     {
       ShapeQualityMetric *met = new VertexConditionNumberQualityMetric;
       MsqError err;
         //Test the metric on a single elemnt patch
       PatchData p1;
       create_one_tri_patch(p1, err);
       MsqVertex* vert1=p1.get_vertex_array(err);
       double first_val=-1;
       bool first_bool=met->evaluate_vertex(p1,&vert1[0],first_val,err);
       if(pF)
         PRINT_INFO("\nVertex Condition No. Metric ideal tri %f", first_val);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(first_val, 1.0, qualTol);
       CPPUNIT_ASSERT(first_bool==true);
       first_bool=met->evaluate_vertex(p1,&vert1[1],first_val,err);
       if(pF)
         PRINT_INFO("\nVertex Condition No. Metric ideal tri %f", first_val);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(first_val, 1.0, qualTol);
       CPPUNIT_ASSERT(first_bool==true);
         //Test on a patch with two ideal tris
       PatchData p2;
       create_two_tri_patch(p2, err);
       MsqVertex* vert2=p2.get_vertex_array(err);
       double second_val;
       bool second_bool=met->evaluate_vertex(p2,&vert2[0],second_val,err);
         //Two neighboring tris with equal area should have local size of 1.0
       if(pF)
         PRINT_INFO("\nVertex Condition No. Metric ideal tris %f", second_val);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(second_val, 1.0, qualTol);
       CPPUNIT_ASSERT(second_bool==true);
     }
  
  void test_untangle_metric()
     {
       UntangleQualityMetric *met = new UntangleBetaQualityMetric(0.0);
       MsqError err;
       double val;
         //Test the metric on a single elemnt patch
       MsqMeshEntity* elem1=idealTri.get_element_array(err);
       bool first_bool=met->evaluate_element(idealTri,&elem1[0],val,err);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(val, 0.0, qualTol);
       CPPUNIT_ASSERT(first_bool==true);
       elem1=idealTet.get_element_array(err);
       first_bool=met->evaluate_element(idealTet,&elem1[0],val,err);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(val, 0.0, qualTol);
       CPPUNIT_ASSERT(first_bool==true);
         //elem1=invertedTri.get_element_array(err);
         //first_bool=met->evaluate_element(invertedTri,&elem1[0],val,err);
         //std::cout<<"\nINVERTED TRI "<<val<<"\n";
         //CPPUNIT_ASSERT_DOUBLES_EQUAL(first_val, 0.0, MSQ_MIN);
         //CPPUNIT_ASSERT(first_bool==true);
       elem1=invertedTet.get_element_array(err);
       first_bool=met->evaluate_element(invertedTet,&elem1[0],val,err);
         //std::cout<<"\nINVERTED TET "<<val<<"\n";
         //Michael:: double check to make sure that 2.0 is correct here.
       CPPUNIT_ASSERT_DOUBLES_EQUAL(val, 2.0, qualTol);
       CPPUNIT_ASSERT(first_bool==true);
     }  
  
  void test_averaging_method()
     {
         //USE QUAD's
     MsqError err;
     double val;
     MsqMeshEntity* elems;
     elems=quadPatch.get_element_array(err);
     MSQ_CHKERR(err);
     ShapeQualityMetric *met = new MeanRatioQualityMetric;
       //Check mean ratio of ideal quad
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
     met->set_averaging_method(QualityMetric::GEOMETRIC, err);
       //Check mean ratio of ideal quad GEOMETRIC
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
     met->set_averaging_method(QualityMetric::HARMONIC, err);
       //Check mean ratio of ideal quad HARMONIC
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
     met->set_averaging_method(QualityMetric::LINEAR, err);
       //Check mean ratio of ideal quad LINEAR
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
     met->set_averaging_method(QualityMetric::MAXIMUM, err);
       //Check mean ratio of ideal quad MAXIMUM
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
     met->set_averaging_method(QualityMetric::MINIMUM, err);
       //Check mean ratio of ideal quad MINIMUM
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
     met->set_averaging_method(QualityMetric::RMS, err);
       //Check mean ratio of ideal quad RMS
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
     met->set_averaging_method(QualityMetric::SUM, err);
       //Check mean ratio of ideal SUM (NOTICE:: should be 4.0)
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,4.0,qualTol);
     met->set_averaging_method(QualityMetric::MAX_OVER_MIN, err);
       //Check mean ratio of ideal MAX_OVER_MIN (NOTICE:: should be 1.0)
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
     
     met->set_averaging_method(QualityMetric::MAX_MINUS_MIN, err);
       //Check mean ratio of ideal MAX_MINUS_MIN (NOTICE:: should be 0.0)
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,0.0,qualTol);
     met->set_averaging_method(QualityMetric::STANDARD_DEVIATION, err);
       //Check mean ratio of ideal STANDARD_DEVIATION (NOTICE:: should be 0.0)
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,0.0,qualTol);
     
     met->set_averaging_method(QualityMetric::SUM_OF_RATIOS_SQUARED, err);
       //Check mean ratio of ideal SUM_OF_RATIOS_SQR (NOTICE:: should be 1.0)
     met->evaluate_element(quadPatch,&elems[0],val,err);MSQ_CHKERR(err);
     CPPUNIT_ASSERT_DOUBLES_EQUAL(val,1.0,qualTol);
   }

  void test_mean_ratio_gradient(PatchData &pd)
  {
    MsqError err;
    Vector3D* grad_num = new Vector3D[2];
    Vector3D* grad_ana = new Vector3D[2];
    double metric_value;
    bool valid;

    MsqMeshEntity* elems = pd.get_element_array(err);MSQ_CHKERR(err);
    MsqVertex* vertices =  pd.get_vertex_array(err);MSQ_CHKERR(err);

    std::vector<size_t> bad_elem_vertex_indices;
    elems[1].get_vertex_indices(bad_elem_vertex_indices);
    MsqVertex* two_vtces[2];
    two_vtces[0] = &vertices[bad_elem_vertex_indices[0]];
    two_vtces[1] = &vertices[bad_elem_vertex_indices[2]];
    
    // creates a mean ratio quality metric ...
    ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
    mean_ratio->set_averaging_method(QualityMetric::SUM, err);

    mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
    valid = mean_ratio->compute_element_gradient (pd, &elems[1], two_vtces,
                                          grad_num, 2, metric_value, err); MSQ_CHKERR(err);
    CPPUNIT_ASSERT(valid);
    
//     std::cout << "NUMERICAL GRADIENT\n";
//     for (int i=0; i<2; ++i)
//        for (int j=0; j<3; ++j)
//          std::cout << grad_num[i][j] << std::endl;

    mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
    valid = mean_ratio->compute_element_gradient (pd, &elems[1], two_vtces,
                                          grad_ana, 2, metric_value, err); MSQ_CHKERR(err);
    CPPUNIT_ASSERT(valid);
//     std::cout << "ANALYTICAL GRADIENT\n";
//     for (int i=0; i<2; ++i)
//       for (int j=0; j<3; ++j)
//         std::cout << grad_ana[i][j] << std::endl;

    for (int i=0; i<2; ++i)
      for (int j=0; j<3; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(grad_num[i][j], grad_ana[i][j], 0.001);
    

    // same test, but free vertices order differ from vertices order in element. 
    two_vtces[0] = &vertices[bad_elem_vertex_indices[2]];
    two_vtces[1] = &vertices[bad_elem_vertex_indices[0]];
    
    mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
    valid = mean_ratio->compute_element_gradient (pd, &elems[1], two_vtces,
                                          grad_num, 2, metric_value, err); MSQ_CHKERR(err);
    CPPUNIT_ASSERT(valid);

    mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
    valid = mean_ratio->compute_element_gradient (pd, &elems[1], two_vtces,
                                          grad_ana, 2, metric_value, err); MSQ_CHKERR(err);
    CPPUNIT_ASSERT(valid);

    for (int i=0; i<2; ++i)
      for (int j=0; j<3; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(grad_num[i][j], grad_ana[i][j], 0.001);

    delete grad_num;
    delete grad_ana;
  }

  void test_mean_ratio_tri_gradient_planar()
  {
    triPatch.vertex_by_index(3).z(0.);
    test_mean_ratio_gradient(triPatch);
  }
       
  void test_mean_ratio_tri_gradient_nonplanar()
  {
    test_mean_ratio_gradient(triPatch);
  }
       
  void test_mean_ratio_quad_gradient_planar()
  {
    quadPatch.vertex_by_index(4).z(0.);
    quadPatch.vertex_by_index(5).z(0.);
    test_mean_ratio_gradient(quadPatch);
  }
       
  void test_mean_ratio_quad_gradient_nonplanar()
  {
    test_mean_ratio_gradient(quadPatch);
  }
       
  void test_mean_ratio_hex_gradient()
  {
    test_mean_ratio_gradient(hexPatch);
  }
       
  void test_mean_ratio_tet_gradient()
  {
    test_mean_ratio_gradient(tetPatch);
  }
  /*! This tests the QualityMetric hessian, comparing analytical
      and numerical versions. 
      
      \param pd: this PatchData must have at least two arguments.
      \param *met: pointer to a metric which will be used in the
      test.  NOTE:  the test may change the current type of gradient
      or hessian evaluations used for the metric.
  */
  void test_metric_hessian(PatchData &pd, QualityMetric *met)
  {
    MsqError err;
    int max_nve = MSQ_MAX_NUM_VERT_PER_ENT;
    Vector3D* grad_num = new Vector3D[max_nve];
    Vector3D* grad_ana = new Vector3D[max_nve];
    Matrix3D* hessian_num = new Matrix3D[max_nve*(max_nve+1)/2];
    Matrix3D* hessian_ana = new Matrix3D[max_nve*(max_nve+1)/2];
    double metric_value;
    double metric_value2;
    

    MsqMeshEntity* elems = pd.get_element_array(err);MSQ_CHKERR(err);
    MsqVertex* vertices =  pd.get_vertex_array(err);MSQ_CHKERR(err);

    std::vector<size_t> elem_vtx_indices;
    elems[1].get_vertex_indices(elem_vtx_indices);
    int nve = elem_vtx_indices.size(); // number of vertices in element.
    MsqVertex** all_vtces = new MsqVertex*[nve];
    for (int i=0; i<nve; ++i) {
      all_vtces[i] = &vertices[elem_vtx_indices[i]];
    }

    all_vtces[0] = &vertices[elem_vtx_indices[0]];
    all_vtces[1] = &vertices[elem_vtx_indices[2]];
    all_vtces[2] = NULL;
    met->set_hessian_type(QualityMetric::NUMERICAL_HESSIAN);
    bool ret_bool=met->compute_element_hessian(pd, &elems[1],
                                               all_vtces, grad_num,
                                               hessian_num, 2,
                                               metric_value,err);
    MSQ_CHKERR(err);
    CPPUNIT_ASSERT(ret_bool==true);
//     std::cout << "GRADIENT for element with two free vertices.\n";
//     for (int i=0; i<4; ++i)
//       for (int j=0; j<3; ++j)
//         std::cout << grad_num[i][j] << std::endl;

//     std::cout << "NUMERICAL HESSIAN for element with two free vertices.\n";
//     for (int i=0; i<nve*(nve+1)/2; ++i)
//          std::cout << hessian_num[i] << std::endl;

    met->set_hessian_type(QualityMetric::ANALYTICAL_HESSIAN);
    ret_bool=met->compute_element_hessian(pd, &elems[1], all_vtces,
                                          grad_ana, hessian_ana, 2,
                                          metric_value2,err);
    MSQ_CHKERR(err);
    CPPUNIT_ASSERT(ret_bool==true);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(metric_value, metric_value2, 0.001); 
//     std::cout << "ANALYTICAL HESSIAN for element with two free vertices.\n";
//     for (int i=0; i<nve*(nve+1)/2; ++i)
//         std::cout << hessian_ana[i] << std::endl;

    // test returned gradients
    for (int m=0; m<nve; ++m)
      for (int i=0; i<3; ++i)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(grad_num[m][i], grad_ana[m][i], 0.001);    

    // test returned Hessians
    for (int m=0; m<nve*(nve+1)/2; ++m)
      for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(hessian_num[m][i][j], hessian_ana[m][i][j], 0.02);


    
    delete[] all_vtces;
    delete[] grad_num;
    delete[] grad_ana;
    delete[] hessian_num;
    delete[] hessian_ana;
  }
  

  
  /*! This tests the QualityMetric hessian, comparing analytical
      and numerical versions. Two comparisons are performed, one for
      elements with free vertices only, and one for an element that
      includes fixed vertices.
      
      \param pd: this PatchData must have at least two arguments.
  */
  void test_mean_ratio_hessian(PatchData &pd)
  {
    MsqError err;
    int max_nve = MSQ_MAX_NUM_VERT_PER_ENT;
    Vector3D* grad_num = new Vector3D[max_nve];
    Vector3D* grad_ana = new Vector3D[max_nve];
    Matrix3D* hessian_num = new Matrix3D[max_nve*(max_nve+1)/2];
    Matrix3D* hessian_ana = new Matrix3D[max_nve*(max_nve+1)/2];
    double metric_value;

    MsqMeshEntity* elems = pd.get_element_array(err);MSQ_CHKERR(err);
    MsqVertex* vertices =  pd.get_vertex_array(err);MSQ_CHKERR(err);

    std::vector<size_t> elem_vtx_indices;
    elems[1].get_vertex_indices(elem_vtx_indices);
    int nve = elem_vtx_indices.size(); // number of vertices in element.
    MsqVertex** all_vtces = new MsqVertex*[nve];
    for (int i=0; i<nve; ++i) {
      all_vtces[i] = &vertices[elem_vtx_indices[i]];
    }

    // 1 **** test with all vertices free
    // creates a mean ratio quality metric ...
    ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
    

//    mean_ratio->set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
    mean_ratio->set_hessian_type(QualityMetric::NUMERICAL_HESSIAN);
    mean_ratio->compute_element_hessian(pd, &elems[1], all_vtces,
                                        grad_num, hessian_num, nve, metric_value,
                                        err); MSQ_CHKERR(err);

//     std::cout << "GRADIENT for element with all  vertices free.\n";
//     for (int i=0; i<nve; ++i)
//       for (int j=0; j<3; ++j)
//         std::cout << grad_num[i][j] << std::endl;

//     std::cout << "NUMERICAL HESSIAN for element with all  vertices free.\n";
//     for (int i=0; i<nve*(nve+1)/2; ++i)
//          std::cout << hessian_num[i] << std::endl;

    mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
    mean_ratio->set_hessian_type(QualityMetric::ANALYTICAL_HESSIAN);
    mean_ratio->compute_element_hessian(pd, &elems[1], all_vtces,
                                        grad_ana, hessian_ana, nve, metric_value,
                                        err); MSQ_CHKERR(err);

//     std::cout << "ANALYTICAL HESSIAN for element with all  vertices free.\n";
//     for (int i=0; i<nve*(nve+1)/2; ++i)
//         std::cout << hessian_ana[i] << std::endl;

    // test returned gradients
    for (int m=0; m<nve; ++m)
      for (int i=0; i<3; ++i)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(grad_num[m][i], grad_ana[m][i], 0.001);    

    // test returned Hessians
    for (int m=0; m<nve*(nve+1)/2; ++m)
      for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j){
            //PRINT_INFO("\nm=%i,i=%i,j=%i",m,i,j);
            //PRINT_INFO("\nNumerical = %f, Analytical = %f",hessian_num[m][i][j], hessian_ana[m][i][j]);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(hessian_num[m][i][j], hessian_ana[m][i][j], 0.003);
        }
    
    
    // 2 **** same test as 1, but gives the free vertices in an order
    //        different than the order within the elements.
    //        Test check that an error is set. 

    // swaps free vertices 0 and 2.
    MsqVertex* swap;
    swap = all_vtces[2];
    all_vtces[2] = all_vtces[0];
    all_vtces[0] = swap;
    
    mean_ratio->set_hessian_type(QualityMetric::ANALYTICAL_HESSIAN);
    bool ret_bool = mean_ratio->compute_element_hessian(pd, &elems[1],
                                                        all_vtces,
                                                        grad_ana, hessian_ana,
                                                        nve, metric_value,
                                                        err);
    
    CPPUNIT_ASSERT(ret_bool == false);
    CPPUNIT_ASSERT(err.errorOn == true);
    err.reset();

    
    // 3 **** same test as 1, but with only 2 free vertices in the element. 
    all_vtces[0] = &vertices[elem_vtx_indices[0]];
    all_vtces[1] = &vertices[elem_vtx_indices[2]];
    all_vtces[2] = NULL;
    mean_ratio->set_hessian_type(QualityMetric::NUMERICAL_HESSIAN);
    ret_bool=mean_ratio->compute_element_hessian(pd, &elems[1], all_vtces,
                                                 grad_num, hessian_num, 2,
                                                 metric_value,
                                                 err); MSQ_CHKERR(err);
    CPPUNIT_ASSERT(ret_bool==true);
//     std::cout << "GRADIENT for element with two free vertices.\n";
//     for (int i=0; i<4; ++i)
//       for (int j=0; j<3; ++j)
//         std::cout << grad_num[i][j] << std::endl;

//     std::cout << "NUMERICAL HESSIAN for element with two free vertices.\n";
//     for (int i=0; i<nve*(nve+1)/2; ++i)
//          std::cout << hessian_num[i] << std::endl;

    mean_ratio->set_hessian_type(QualityMetric::ANALYTICAL_HESSIAN);
    ret_bool=mean_ratio->compute_element_hessian(pd, &elems[1], all_vtces,
                                                 grad_ana, hessian_ana, 2,
                                                 metric_value,
                                                 err); MSQ_CHKERR(err);
    CPPUNIT_ASSERT(ret_bool==true);
//     std::cout << "ANALYTICAL HESSIAN for element with two free vertices.\n";
//     for (int i=0; i<nve*(nve+1)/2; ++i)
//         std::cout << hessian_ana[i] << std::endl;

    // test returned gradients
    for (int m=0; m<nve; ++m)
      for (int i=0; i<3; ++i)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(grad_num[m][i], grad_ana[m][i], 0.001);    

    // test returned Hessians
    for (int m=0; m<nve*(nve+1)/2; ++m)
      for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(hessian_num[m][i][j], hessian_ana[m][i][j], 0.003);

       
 
    delete[] all_vtces;
    delete[] grad_num;
    delete[] grad_ana;
    delete[] hessian_num;
    delete[] hessian_ana;
  }
       
  void test_mean_ratio_tri_hessian()
    {
      MsqError err;
      test_mean_ratio_hessian(triPatch);
      QualityMetric* mean_rat = new MeanRatioQualityMetric;
        //make sure the code handles this case correctly
      mean_rat->set_averaging_method(QualityMetric::SUM_SQUARED,
                                     err);
      test_metric_hessian(triPatch, mean_rat);
    }
  
        
  
  void test_mean_ratio_quad_hessian()
  {
    test_mean_ratio_hessian(quadPatch);
  }
  void test_mean_ratio_quad_hessian_linear()
    {
      MsqError err;
      QualityMetric* mean_rat = new MeanRatioQualityMetric;
      mean_rat->set_averaging_method(QualityMetric::LINEAR,
                                   err);
      test_metric_hessian(quadPatch, mean_rat);
      delete mean_rat;
    }
  void test_mean_ratio_quad_hessian_sum_squared()
    {
      MsqError err;
      QualityMetric* mean_rat = new MeanRatioQualityMetric;
      mean_rat->set_averaging_method(QualityMetric::SUM_SQUARED,
                                     err);
      test_metric_hessian(quadPatch, mean_rat);
      delete mean_rat;
    }
   void test_mean_ratio_quad_hessian_rms()
    {
      MsqError err;
      QualityMetric* mean_rat = new MeanRatioQualityMetric;
      mean_rat->set_averaging_method(QualityMetric::RMS,
                                     err);
      test_metric_hessian(quadPatch, mean_rat);
      delete mean_rat;
    }
  void test_mean_ratio_quad_hessian_harmonic()
    {
      MsqError err;
      QualityMetric* mean_rat = new MeanRatioQualityMetric;
      mean_rat->set_averaging_method(QualityMetric::HARMONIC,
                                     err);
      test_metric_hessian(quadPatch, mean_rat);
      delete mean_rat;
    }
  void test_mean_ratio_quad_hessian_hms()
    {
      MsqError err;
      QualityMetric* mean_rat = new MeanRatioQualityMetric;
      mean_rat->set_averaging_method(QualityMetric::HMS,
                                     err);
      test_metric_hessian(quadPatch, mean_rat);
      delete mean_rat;
    }
  void test_mean_ratio_tet_hessian()
  {
    test_mean_ratio_hessian(tetPatch);
  }
  
  void test_mean_ratio_hex_hessian()
  {
    test_mean_ratio_hessian(hexPatch);
     MsqError err;
    QualityMetric* mean_rat = new MeanRatioQualityMetric;
    mean_rat->set_averaging_method(QualityMetric::SUM_SQUARED,
                                   err);
    test_metric_hessian(hexPatch, mean_rat);
  }
 void test_mean_ratio_hex_hessian_linear()
    {
      MsqError err;
      QualityMetric* mean_rat = new MeanRatioQualityMetric;
      mean_rat->set_averaging_method(QualityMetric::LINEAR,
                                     err);
      test_metric_hessian(hexPatch, mean_rat);
      delete mean_rat;
    }  
 void test_mean_ratio_hex_hessian_sum_squared()
    {
      MsqError err;
      QualityMetric* mean_rat = new MeanRatioQualityMetric;
      mean_rat->set_averaging_method(QualityMetric::SUM_SQUARED,
                                     err);
      test_metric_hessian(hexPatch, mean_rat);
      delete mean_rat;
    }  
  void test_mean_ratio_hex_hessian_rms()
    {
      MsqError err;
      QualityMetric* mean_rat = new MeanRatioQualityMetric;
      mean_rat->set_averaging_method(QualityMetric::RMS,
                                     err);
      test_metric_hessian(hexPatch, mean_rat);
      delete mean_rat;
    }
 void test_mean_ratio_hex_hessian_harmonic()
    {
      MsqError err;
      QualityMetric* mean_rat = new MeanRatioQualityMetric;
      mean_rat->set_averaging_method(QualityMetric::HARMONIC,
                                     err);
      test_metric_hessian(hexPatch, mean_rat);
      delete mean_rat;
    }
  void test_mean_ratio_hex_hessian_hms()
    {
      MsqError err;
      QualityMetric* mean_rat = new MeanRatioQualityMetric;
      mean_rat->set_averaging_method(QualityMetric::HMS,
                                     err);
      test_metric_hessian(hexPatch, mean_rat);
      delete mean_rat;
    }
    /*!       
      \param pd: this PatchData must have at least two elements.
  */
  void test_mean_ratio_grad_from_hessian(PatchData &pd)
  {
    MsqError err;
    int max_nve = MSQ_MAX_NUM_VERT_PER_ENT;
    Vector3D* grad1 = new Vector3D[max_nve];
    Vector3D* grad2 = new Vector3D[max_nve];
    Matrix3D* hessian = new Matrix3D[max_nve*(max_nve+1)/2];
    double QM_val1, QM_val2;
    bool valid;
    
    MsqMeshEntity* elems = pd.get_element_array(err);MSQ_CHKERR(err);
    MsqVertex* vertices =  pd.get_vertex_array(err);MSQ_CHKERR(err);

    std::vector<size_t> elem_vtx_indices;
    elems[1].get_vertex_indices(elem_vtx_indices);
    int nve = elem_vtx_indices.size(); // number of vertices in element.
    MsqVertex** all_vtces = new MsqVertex*[nve];
    for (int i=0; i<nve; ++i) {
      all_vtces[i] = &vertices[elem_vtx_indices[i]];
    }

    // 1 **** test with all vertices free
    // creates a mean ratio quality metric ...
    ShapeQualityMetric* mean_ratio = new MeanRatioQualityMetric;
    mean_ratio->set_averaging_method(QualityMetric::SUM, err); MSQ_CHKERR(err);
    mean_ratio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
    mean_ratio->set_hessian_type(QualityMetric::ANALYTICAL_HESSIAN);


    valid = mean_ratio->compute_element_gradient (pd, &elems[1],
                                                  all_vtces, grad1,
                                                  nve, QM_val1, err);
    MSQ_CHKERR(err); CPPUNIT_ASSERT(valid);

    
    valid = mean_ratio->compute_element_hessian(pd, &elems[1], all_vtces,
                                        grad2, hessian,
                                        nve, QM_val2,
                                        err); MSQ_CHKERR(err);
    CPPUNIT_ASSERT(valid);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(QM_val1, QM_val2, 1e-12);

//     cout << "Gradient from compute_gradient()\n";
//     for (int i=0; i<nve; ++i)
//       cout << grad1[i];

//     cout << "\nGradient from compute_hessian()\n";
//     for (int i=0; i<nve; ++i)
//       cout << grad2[i];

    
    // test gradients
    for (int i=0; i<nve; ++i)
      for (int j=0; j<3; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(grad1[i][j], grad2[i][j], 1e-12);
  }

  void test_mean_ratio_tri_grad_from_hessian()
  {
    test_mean_ratio_grad_from_hessian(triPatch);
  }
  
  void test_mean_ratio_quad_grad_from_hessian()
  {
    test_mean_ratio_grad_from_hessian(quadPatch);
  }
  
  void test_mean_ratio_tet_grad_from_hessian()
  {
    test_mean_ratio_grad_from_hessian(tetPatch);
  }
  
  void test_mean_ratio_hex_grad_from_hessian()
  {
    test_mean_ratio_grad_from_hessian(hexPatch);
  }
  
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QualityMetricTest, "QualityMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(QualityMetricTest, "Unit");
