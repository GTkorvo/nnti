// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 12-Nov-02 at 18:05:56
//  LAST-MOD:  7-May-03 at 13:28:47 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MsqMeshEntityTest.cpp

Unit testing of various functions in the MsqMeshEntity class.

\author Michael Brewer
\author Thomas Leurent

 */
// DESCRIP-END.
//


#include "MsqMeshEntity.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"
#include "PatchDataInstances.hpp"
#include <math.h>
#include "MsqMessage.hpp"
#include <iostream>
#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

using namespace Mesquite;
using std::cout;
using std::endl;

class MsqMeshEntityTest : public CppUnit::TestFixture
{

private:
  CPPUNIT_TEST_SUITE(MsqMeshEntityTest);
  CPPUNIT_TEST (test_hex_vertices);
  CPPUNIT_TEST (test_centroid_tri);
  CPPUNIT_TEST (test_centroid_quad);
  CPPUNIT_TEST (test_centroid_hex);
  CPPUNIT_TEST (test_compute_weigted_jacobian_ideal_hex);
  CPPUNIT_TEST (test_compute_weigted_jacobian_ideal_tet);
  CPPUNIT_TEST (test_compute_weigted_jacobian_ideal_quad);
  CPPUNIT_TEST (test_compute_weigted_jacobian_ideal_tri);
  CPPUNIT_TEST (test_unsigned_area);
  CPPUNIT_TEST_SUITE_END();

private:
  PatchData oneHexPatch;
  PatchData oneTetPatch;
  PatchData oneQuadPatch;
  PatchData oneTriPatch;
  Vector3D e1, e2, e3;
  double tolEps;
public:
  void setUp()
  {
    tolEps=1.e-12;
    
    // sets up the unit vectors
    e1.set(1,0,0);
    e2.set(0,1,0);
    e3.set(0,0,1);

    MsqError err;

    // creates empty Patch
    create_one_hex_patch(oneHexPatch, err); MSQ_CHKERR(err);
    create_one_tet_patch(oneTetPatch, err); MSQ_CHKERR(err);
    create_one_tri_patch(oneTriPatch, err); MSQ_CHKERR(err);
    create_one_quad_patch(oneQuadPatch, err); MSQ_CHKERR(err);
  }

  void tearDown()
  {

  }
  
public:
  MsqMeshEntityTest()
  {  }
  
  void test_hex_vertices()
  {
    MsqError err;
    // prints out the vertices.
    MsqVertex* ideal_vertices = oneHexPatch.get_vertex_array(err); MSQ_CHKERR(err);
    size_t num_vtx = oneHexPatch.num_vertices();
    CPPUNIT_ASSERT_EQUAL(size_t(8), num_vtx);
    
    MsqVertex vtx;

    vtx.set(1,1,1);
    CPPUNIT_ASSERT_EQUAL(vtx, ideal_vertices[0]);
    
    vtx.set(2,2,2);
    CPPUNIT_ASSERT_EQUAL(vtx, ideal_vertices[6]);
    
    vtx.set(1,2,2);
    CPPUNIT_ASSERT_EQUAL(vtx, ideal_vertices[7]);
  }

  //! test the centroid of the first element in the Patch
  void test_centroid(PatchData &pd, Vector3D &correct)
  {
    MsqError err;
    double eps = 1e-6;
    Vector3D centroid;

    MsqMeshEntity* elem = pd.get_element_array(err); MSQ_CHKERR(err);
    elem->get_centroid(centroid, pd, err); MSQ_CHKERR(err);

//     cout << "centroid: "<< centroid <<endl; 
//     cout << "correct: "<< correct <<endl; 

    for (int i=0; i<3; ++i)
      CPPUNIT_ASSERT_DOUBLES_EQUAL( centroid[i], correct[i], eps);
  }

  void test_centroid_tri()
  {
    Vector3D correct(1.5, 1+1/(2*sqrt(3)), 1.0);    
    test_centroid(oneTriPatch, correct);
  }

  void test_centroid_quad()
  {
    Vector3D correct(1.5, 1.5, 1.0);    
    test_centroid(oneQuadPatch, correct);
  }

  void test_centroid_hex()
  {
    Vector3D correct(1.5, 1.5, 1.5);    
    test_centroid(oneHexPatch, correct);
  }

  void test_compute_weigted_jacobian_ideal_hex()
  {
    MsqError err;
    //determinate of jacobs
    double deter;
    MsqMeshEntity* hex = oneHexPatch.get_element_array(err); MSQ_CHKERR(err);
    // get the ideal tets sample points 
    std::vector<Vector3D> sample_points;
    hex->get_sample_points(QualityMetric::ELEMENT_VERTICES, sample_points, err); MSQ_CHKERR(err);
      // and get the jacobian vectors (should be the identity matrix).
      //for the first sample point (0,0,0)
    short num_jac_vec;
    std::vector<Vector3D>::iterator sp;
    sp=sample_points.begin();
    Vector3D jacobian_vectors[3];
    hex->compute_weighted_jacobian(oneHexPatch, *sp, jacobian_vectors,
                                   num_jac_vec , err); MSQ_CHKERR(err);
      //Make sure we have the identity
    CPPUNIT_ASSERT((jacobian_vectors[0]-e1).length()<tolEps);
    CPPUNIT_ASSERT((jacobian_vectors[1]-e2).length()<tolEps);
    CPPUNIT_ASSERT((jacobian_vectors[2]-e3).length()<tolEps);
    int counter=0;
    while(sp!=sample_points.end()){
        //get next sample point
      hex->compute_weighted_jacobian(oneHexPatch, *sp, jacobian_vectors,
                                     num_jac_vec , err); MSQ_CHKERR(err);
        //make sure that our jacobian has determinat one
      deter=jacobian_vectors[0]%(jacobian_vectors[1]*jacobian_vectors[2]);
        //<<"\n\nInside test tet c_weighted_j deter = "<< deter <<"\n";
      CPPUNIT_ASSERT(fabs(deter-1.0)<tolEps);
      ++sp;
      ++counter;
    }
      //make sure that we had eight sample points.
    CPPUNIT_ASSERT(counter==8);
  }

       
  void test_compute_weigted_jacobian_ideal_tet()
  {
    MsqError err;
      //determinate of jacobs
    double deter;
    MsqMeshEntity* tet = oneTetPatch.get_element_array(err); MSQ_CHKERR(err);
      // get the ideal tets sample points 
    std::vector<Vector3D> sample_points;
    tet->get_sample_points(QualityMetric::ELEMENT_VERTICES, sample_points, err); MSQ_CHKERR(err);
    // and get the jacobian vectors (should be the identity matrix).
    short num_jac_vec;
    std::vector<Vector3D>::iterator sp;
    sp=sample_points.begin();
    Vector3D jacobian_vectors[3];
    tet->compute_weighted_jacobian(oneTetPatch, *sp, jacobian_vectors,
                                   num_jac_vec , err); MSQ_CHKERR(err);

    CPPUNIT_ASSERT((jacobian_vectors[0]-e1).length()<tolEps);
    CPPUNIT_ASSERT((jacobian_vectors[1]-e2).length()<tolEps);
    CPPUNIT_ASSERT((jacobian_vectors[2]-e3).length()<tolEps);
      //loop over sample points
    int counter=0;
    while(sp!=sample_points.end()){
        //get next sample point
      
      tet->compute_weighted_jacobian(oneTetPatch, *sp, jacobian_vectors,
                                     num_jac_vec , err); MSQ_CHKERR(err);
      deter=jacobian_vectors[0]%(jacobian_vectors[1]*jacobian_vectors[2]);
        //<<"\n\nInside test tet c_weighted_j deter = "<< deter <<"\n";
      CPPUNIT_ASSERT(fabs(deter-1.0)<tolEps);
      ++sp;
      ++counter;
    }
    CPPUNIT_ASSERT(counter==4);
    
  }
void test_compute_weigted_jacobian_ideal_quad()
  {
    MsqError err;
      //determinate of jacobs
    double deter;
    MsqMeshEntity* quad = oneQuadPatch.get_element_array(err); MSQ_CHKERR(err);
    // get the ideal quads sample points 
    std::vector<Vector3D> sample_points;
    quad->get_sample_points(QualityMetric::ELEMENT_VERTICES, sample_points, err); MSQ_CHKERR(err);
    // and get the jacobian vectors (should be the identity matrix).
    short num_jac_vec;
    std::vector<Vector3D>::iterator sp;
    sp=sample_points.begin();
    Vector3D jacobian_vectors[3];
    quad->compute_weighted_jacobian(oneQuadPatch, *sp, jacobian_vectors,
                                   num_jac_vec , err); MSQ_CHKERR(err);
    CPPUNIT_ASSERT((jacobian_vectors[0]-e1).length()<tolEps);
    CPPUNIT_ASSERT((jacobian_vectors[1]-e2).length()<tolEps);
      //loop over sample points
    int counter=0;
    while(sp!=sample_points.end()){
        //get next sample point
      
      quad->compute_weighted_jacobian(oneQuadPatch, *sp, jacobian_vectors,
                                     num_jac_vec , err); MSQ_CHKERR(err);
        //deter is not the determinant in quad case
      deter=(jacobian_vectors[0]*jacobian_vectors[1]).length();
      CPPUNIT_ASSERT(fabs(deter-1.0)<tolEps);
      ++sp;
      ++counter;
    }
    CPPUNIT_ASSERT(counter==4);
    
  }
  
void test_compute_weigted_jacobian_ideal_tri()
  {
    MsqError err;
      //determinate of jacobs
    double deter;
    MsqMeshEntity* tri = oneTriPatch.get_element_array(err); MSQ_CHKERR(err);
    // get the ideal tris sample points 
    std::vector<Vector3D> sample_points;
    tri->get_sample_points(QualityMetric::ELEMENT_VERTICES, sample_points, err); MSQ_CHKERR(err);
    // and get the jacobian vectors (should be the identity matrix).
    short num_jac_vec;
    std::vector<Vector3D>::iterator sp;
    sp=sample_points.begin();
    Vector3D jacobian_vectors[3];
    tri->compute_weighted_jacobian(oneTriPatch, *sp, jacobian_vectors,
                                   num_jac_vec , err); MSQ_CHKERR(err);
    CPPUNIT_ASSERT((jacobian_vectors[0]-e1).length()<tolEps);
    CPPUNIT_ASSERT((jacobian_vectors[1]-e2).length()<tolEps);
      //loop over sample points
    int counter=0;
    while(sp!=sample_points.end()){
        //get next sample point
      
      tri->compute_weighted_jacobian(oneTriPatch, *sp, jacobian_vectors,
                                     num_jac_vec , err); MSQ_CHKERR(err);
      deter=(jacobian_vectors[0]*jacobian_vectors[1]).length();
      CPPUNIT_ASSERT(fabs(deter-1.0)<tolEps);
      ++sp;
      ++counter;
    }
    CPPUNIT_ASSERT(counter==3);
    
  }

  void test_unsigned_area()
     {
       MsqError err;
       MsqMeshEntity* tri = oneTriPatch.get_element_array(err);
       MSQ_CHKERR(err);
       CPPUNIT_ASSERT(fabs(tri->compute_unsigned_area(oneTriPatch,err)
                           -(sqrt(3.0)/4.0)) < tolEps);
       MsqMeshEntity* quad = oneQuadPatch.get_element_array(err);
       MSQ_CHKERR(err);
       CPPUNIT_ASSERT(fabs(quad->compute_unsigned_area(oneQuadPatch,err)
                           -1.0) < tolEps);
       MsqMeshEntity* hex = oneHexPatch.get_element_array(err);
       MSQ_CHKERR(err);
         //PRINT_INFO("\nHEX _ _ _ %f\n",hex->compute_unsigned_volume(oneHexPatch,err));
       CPPUNIT_ASSERT(fabs(hex->compute_unsigned_volume(oneHexPatch,err)
                           -1.0) < tolEps);
     }
  
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqMeshEntityTest, "MsqMeshEntityTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqMeshEntityTest, "Unit");
