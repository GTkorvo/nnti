// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
/*!
  \file   MeshInterface.hpp
  \brief  

  \author Thomas Leurent
  \date   2003-04-17
*/

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "MeshSet.hpp"
#include "PatchData.hpp"
#include "MeshImpl.hpp"
#include "Vector3D.hpp"

#include "cppunit/extensions/HelperMacros.h"
#include "cppunit/SignalException.h"

#include <iostream>


using std::cout;
using std::cerr;
using std::endl;

class MeshInterfaceTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(MeshInterfaceTest);
  CPPUNIT_TEST (test_get_geometric_dimension);
  CPPUNIT_TEST (test_vertices);
  CPPUNIT_TEST (test_vertex_is_on_boundary);
  CPPUNIT_WORK_IN_PROGRESS (test_vertex_is_fixed);
  CPPUNIT_TEST (test_vertex_byte);
  CPPUNIT_TEST (test_vertex_get_attached_elements);
  CPPUNIT_TEST (test_elements);
  CPPUNIT_TEST (test_elements_get_attached_vertices);
  CPPUNIT_TEST (test_element_get_topology);
  CPPUNIT_TEST (test_elements_get_topology);
  CPPUNIT_TEST(test_element_get_attached_vertex_indices);
  CPPUNIT_TEST_SUITE_END();

private:
  Mesquite::MeshImpl *mMesh; 
  Mesquite::Mesh::VertexHandle* mVertices;
  Mesquite::Mesh::ElementHandle* mElements;
  size_t nbVert;
  size_t nbElem;
  Mesquite::MsqError mErr;
  
public:
   /* Automatically called by CppUnit before each test function. */
  void setUp()
  {
      // Read a VTK file -- 1 triangle flanked by 1 quad on each side (1 tri + 3 quads)
    mMesh = new Mesquite::MeshImpl;
    mMesh->read_vtk("../../meshFiles/2D/VTK/hybrid_3quad_1tri.vtk", mErr);
    MSQ_CHKERR(mErr);

    // Gets an array of vertices handles
    nbVert = mMesh->get_total_vertex_count(mErr); MSQ_CHKERR(mErr);
    mVertices = new Mesquite::Mesh::VertexHandle[nbVert];
    mMesh->get_all_vertices(mVertices, nbVert, mErr); MSQ_CHKERR(mErr);

    // Gets an array of element handles
    nbElem = mMesh->get_total_element_count(mErr); MSQ_CHKERR(mErr);
    mElements = new Mesquite::Mesh::ElementHandle[nbElem];
    mMesh->get_all_elements(mElements, nbElem, mErr); MSQ_CHKERR(mErr);

  }
  
    // Automatically called by CppUnit after each test function.
  void tearDown()
  {

    mMesh->release();
    delete [] mVertices;
    delete [] mElements;
    MSQ_CHKERR(mErr);
    mErr.reset();
  }
  
public:
  MeshInterfaceTest()
    {}
  
  void test_get_geometric_dimension()
  {
    int d = mMesh->get_geometric_dimension(mErr);
    CPPUNIT_ASSERT_EQUAL(d,3);
  }

  
#undef __FUNC__
#define __FUNC__ "MeshInterfaceTest::test_vertices"
  void test_vertices()
  {
    CPPUNIT_ASSERT_EQUAL(9,(int)nbVert);

    Mesquite::Vector3D coords;
    Mesquite::Vector3D* correct_coords = new Mesquite::Vector3D[nbVert];
    correct_coords[0].set(1,0,0);
    correct_coords[1].set(0,1.732,0);
    correct_coords[2].set(-1,0,0);
    correct_coords[3].set(-1,-2,0);
    correct_coords[4].set(1,-2,0);
    correct_coords[5].set(2.732,1,0);
    correct_coords[6].set(1.732,2.732,0);
    correct_coords[7].set(-1.732,2.732,0);
    correct_coords[8].set(-2.732,1,0);

    for (size_t i=0; i<nbVert; ++i) {
      mMesh->vertex_get_coordinates(mVertices[i], coords, mErr); MSQ_CHKERR(mErr);
      for (int j=0; j<3; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(coords[j], correct_coords[i][j], .01);
    }

    coords.set(2.,3.,4.);
    mMesh->vertex_set_coordinates(mVertices[3], coords, mErr); MSQ_CHKERR(mErr);
    Mesquite::Vector3D coords_2;
    mMesh->vertex_get_coordinates(mVertices[3], coords_2, mErr); MSQ_CHKERR(mErr);
    for (int j=0; j<3; ++j)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(coords[j], coords_2[j], 1e-6);
    
    delete [] correct_coords;
  }

#undef __FUNC__
#define __FUNC__ "MeshInterfaceTest::test_vertex_is_on_boundary"
  void test_vertex_is_on_boundary()
  {
    bool correct_boundary[9] = {false, false, false, true, true, true, true, true, true};
    for (size_t i=0; i<nbVert; ++i) {
      bool on_bnd = mMesh->vertex_is_on_boundary(mVertices[i], mErr); MSQ_CHKERR(mErr);
      CPPUNIT_ASSERT(on_bnd == correct_boundary[i]);
    }
  }

#undef __FUNC__
#define __FUNC__ "MeshInterfaceTest::test_vertex_is_fixed"
  void test_vertex_is_fixed()
  {
    bool correct_fixed[9] = {false, false, false, true, true, true, true, true, true};
    for (size_t i=0; i<nbVert; ++i) {
      bool fixed = mMesh->vertex_is_fixed(mVertices[i], mErr); MSQ_CHKERR(mErr);
//      cout << "fixed["<<i<<"] : "<< fixed << endl;
      CPPUNIT_ASSERT(fixed == correct_fixed[i]);
    }
  }

#undef __FUNC__
#define __FUNC__ "MeshInterfaceTest::test_vertex_byte"
  void test_vertex_byte()
  {
    unsigned char* bytes = new unsigned char[nbVert];
    mMesh->vertices_get_byte(mVertices, bytes, nbVert, mErr); MSQ_CHKERR(mErr);

    // Asserts all vertex bytes are initialised to 0. 
    for (size_t i=0; i<nbVert; ++i)
      CPPUNIT_ASSERT(bytes[i] == 0);

    // Test various vertex byte read / write routines.
    bytes[3] |= 4;
    mMesh->vertices_set_byte(mVertices, bytes, nbVert, mErr); MSQ_CHKERR(mErr);
    mMesh->vertex_set_byte(mVertices[5], 8, mErr); MSQ_CHKERR(mErr);
    unsigned char byte;
    mMesh->vertex_get_byte(mVertices[3], &byte, mErr); MSQ_CHKERR(mErr);
    CPPUNIT_ASSERT(byte == 4);
    mMesh->vertices_get_byte(mVertices, bytes, nbVert, mErr); MSQ_CHKERR(mErr);
    for (size_t i=0; i<nbVert; ++i) {
      if (i==3)
        CPPUNIT_ASSERT(bytes[i] == 4);
      else if (i==5)
        CPPUNIT_ASSERT(bytes[i] == 8);
      else
        CPPUNIT_ASSERT(bytes[i] == 0);
    }

    delete [] bytes;
  }

  
#undef __FUNC__
#define __FUNC__ "MeshInterfaceTest::test_vertex_get_attached_elements"
  void test_vertex_get_attached_elements()
  {

    // checks we have 6 vertices contained in 1 element only
    // and 3 vertices contained in 3 elements. 
    int n1=0;
    int n3=0;
    for (size_t i=0; i<nbVert; ++i) {
      size_t nev = mMesh->vertex_get_attached_element_count(mVertices[i], mErr);
      MSQ_CHKERR(mErr);
      if (nev==1)
        ++n1;
      else if (nev==3)
        ++n3;
      else // failure. 
        CPPUNIT_ASSERT(false);
    }
    CPPUNIT_ASSERT(n1==6);
    CPPUNIT_ASSERT(n3==3);

    // gets the index of a vertex in a corner
    int one_corner_vertex_index=0;
    size_t i=nbVert-1;
    while (!one_corner_vertex_index) {
      size_t nev = mMesh->vertex_get_attached_element_count(mVertices[i], mErr);
      if (nev==1)
        one_corner_vertex_index=i;
      --i;
    }

    // retrieve the attached element.
    Mesquite::Mesh::ElementHandle elem=0;
    mMesh->vertex_get_attached_elements(mVertices[2], &elem, 1, mErr);
    MSQ_CHKERR(mErr);
    CPPUNIT_ASSERT(elem!=0);

  }
  
#undef __FUNC__
#define __FUNC__ "MeshInterfaceTest::test_elements"
  void test_elements()
  {
    CPPUNIT_ASSERT_EQUAL(4,(int)nbElem);

    
  }


#undef __FUNC__
#define __FUNC__ "MeshInterfaceTest::test_elements_get_attached_vertices"
  void test_elements_get_attached_vertices()
  {
    // checks we have 3 elements containing 4 vertices 
    // and 1 element containing 3 vertices. 
    int n3=0;
    int n4=0;
    for (size_t i=0; i<nbElem; ++i) {
      size_t nve = mMesh->element_get_attached_vertex_count(
                                         mElements[i], mErr);
      MSQ_CHKERR(mErr);
      if (nve==3)
        ++n3;
      else if (nve==4)
        ++n4;
      else // failure. 
        CPPUNIT_ASSERT(false);
    }
    CPPUNIT_ASSERT(n3==1); // 1 triangle
    CPPUNIT_ASSERT(n4==3); // 3 quads

    size_t sizeof_vert_handles = nbVert+3; // testing also with extra space handling
    Mesquite::Mesh::VertexHandle* vert_handles =
      new Mesquite::Mesh::VertexHandle[sizeof_vert_handles];
    size_t sizeof_csr_data = 15;
    size_t* csr_data = new size_t[sizeof_csr_data];
    size_t* csr_offsets = new size_t[nbElem+1];
    
    mMesh->elements_get_attached_vertices(mElements, nbElem,
                                          vert_handles, sizeof_vert_handles,
                                          csr_data, sizeof_csr_data,
                                          csr_offsets,
                                          mErr); MSQ_CHKERR(mErr);

//     cout << "nbElem " << nbElem << endl;
//     for (size_t i=0; i<nbElem; ++i) {
//       cout << "mElements["<<i<<"]: " << mElements[i] << endl;
//     }
//     cout << "sizeof_vert_handles " << sizeof_vert_handles << endl;
//     for (size_t i=0; i< sizeof_vert_handles ; ++i) {
//       cout << "vert_handles["<<i<<"]: " << vert_handles[i] << endl;
//     }
//     cout << "sizeof_csr_data " << sizeof_csr_data << endl;
//     for (size_t i=0; i< sizeof_csr_data ; ++i) {
//       cout << "csr_data["<<i<<"]: " << csr_data[i] << endl;
//     }
//     for (size_t i=0; i< nbElem+1 ; ++i) {
//       cout << "csr_offsets["<<i<<"]: " << csr_offsets[i] << endl;
//     }

    // Make sure a handle is returned for each vertex.
    for (int i=0; i<9; ++i) {
      CPPUNIT_ASSERT( vert_handles[i] != 0 );
    }
    
    // Make sure CSR data is valid
    int vtx_repeated_occurence[9] = {0,0,0,0,0,0,0,0,0};
    for (int i=0; i<15; ++i) {
      CPPUNIT_ASSERT( csr_data[i] >=0 && csr_data[i] <= 8 );
      ++vtx_repeated_occurence[csr_data[i]];
    }
    for (int i=0; i<9; ++i) {
      CPPUNIT_ASSERT( vtx_repeated_occurence[i] <= 3 );
    }

    // Makes sure CSR offsets are valid
    CPPUNIT_ASSERT( csr_offsets[0] == 0 );
    CPPUNIT_ASSERT( csr_offsets[1] >=3 && csr_offsets[1]<=12 );
    CPPUNIT_ASSERT( csr_offsets[2] >=3 && csr_offsets[2]<=12 );
    CPPUNIT_ASSERT( csr_offsets[3] >=3 && csr_offsets[3]<=12 );
    CPPUNIT_ASSERT( csr_offsets[4] == 15 );
    
    // When sizeof_csr_data is insufficient, makes sure an error is set 
    sizeof_csr_data = 12;
    mMesh->elements_get_attached_vertices(mElements, nbElem,
                                          vert_handles, sizeof_vert_handles,
                                          csr_data, sizeof_csr_data,
                                          csr_offsets,
                                          mErr);
    CPPUNIT_ASSERT( mErr.errorOn == true );
    mErr.reset();
    delete []vert_handles;
    delete []csr_data;
    delete []csr_offsets;
  }

  
#undef __FUNC__
#define __FUNC__ "MeshInterfaceTest::test_element_get_topology"
  void test_element_get_topology()
  {
    int nb_quads=0;
    int nb_tri=0;
    Mesquite::EntityTopology topo;
    for (size_t i=0; i<nbElem; ++i) {
      topo = mMesh->element_get_topology(mElements[i], mErr);
      MSQ_CHKERR(mErr);
      switch (topo) {
      case Mesquite::TRIANGLE:
        ++nb_tri;
        break;
      case Mesquite::QUADRILATERAL:
        ++nb_quads;
        break;
      default:
        CPPUNIT_FAIL("Topology should be quad or Hex only.");
      }
    }
    CPPUNIT_ASSERT_EQUAL(1,nb_tri);
    CPPUNIT_ASSERT_EQUAL(3,nb_quads);
  }

#undef __FUNC__
#define __FUNC__ "MeshInterfaceTest::test_elements_get_topology"
  void test_elements_get_topology()
  {
    int nb_quads=0;
    int nb_tri=0;
    Mesquite::EntityTopology* topos = new Mesquite::EntityTopology[nbElem];
    mMesh->elements_get_topologies(mElements, topos, nbElem, mErr);
    MSQ_CHKERR(mErr);
    for (size_t i=0; i<nbElem; ++i) {
      switch (topos[i]) {
      case Mesquite::TRIANGLE:
        ++nb_tri;
        break;
      case Mesquite::QUADRILATERAL:
        ++nb_quads;
        break;
      default:
        CPPUNIT_FAIL("Topology should be quad or Hex only.");
      }
    }
    CPPUNIT_ASSERT_EQUAL(1,nb_tri);
    CPPUNIT_ASSERT_EQUAL(3,nb_quads);
    delete []topos;
  }



#undef __FUNC__
#define __FUNC__ "MeshInterfaceTest::test_element_get_attached_vertex_indices"
  void test_element_get_attached_vertex_indices()
  {
    // Find the index of the triangle
    Mesquite::EntityTopology topo=Mesquite::MIXED;
    int tri_index=-1;
    while (topo != Mesquite::TRIANGLE) {
      ++tri_index;
      topo = mMesh->element_get_topology(mElements[tri_index], mErr);
      MSQ_CHKERR(mErr);
    }

    size_t index_array[3];
    mMesh->element_get_attached_vertex_indices(mElements[tri_index],
                                        index_array, 3, mErr);

    // creates list with correct vertices coordinates for the triangle
    std::list<Mesquite::Vector3D> correct_coords;
    correct_coords.push_back(Mesquite::Vector3D(1.,0.,0.));
    correct_coords.push_back(Mesquite::Vector3D(0.,1.732050807,0.));
    correct_coords.push_back(Mesquite::Vector3D(-1.,0.,0.));

    // Creates same list from the mesh implementation
    std::list<Mesquite::Vector3D> tri_coords;
    Mesquite::Vector3D tmp_vec;
    for (size_t i=0; i<3; ++i) {
      mMesh->vertex_get_coordinates(mVertices[i], tmp_vec, mErr);
      tri_coords.push_back(tmp_vec);
    }

    // Makes sure both list contain the same elements (not necessarily in the same order).
    std::list<Mesquite::Vector3D>::iterator correct_iter;
    correct_iter = correct_coords.begin();
    std::list<Mesquite::Vector3D>::iterator tri_iter;
    while (!correct_coords.empty()) {
      //      cout << "correct_iter : " << *correct_iter;
      tri_iter = tri_coords.begin();
      while (tri_iter != tri_coords.end()) {
        //  cout << "tri_iter : " << *tri_iter;
        if (Mesquite::Vector3D::distance_between(*tri_iter, *correct_iter) < 10e-4) {
          //cout << "Erasing received : "<< *tri_iter;
          tri_coords.erase( tri_iter );
          tri_iter = tri_coords.end(); // stops
          //cout << "Erasing correct : "<< *correct_iter;
          correct_coords.erase( correct_iter );
        }
        else
          ++tri_iter;
      }
      correct_iter = correct_coords.begin();
    }
    CPPUNIT_ASSERT(tri_coords.empty());
    CPPUNIT_ASSERT(correct_coords.empty());
  }

};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MeshInterfaceTest, "MeshInterfaceTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MeshInterfaceTest, "Unit");
