/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retian certain rights to this software.

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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file VertexConditionNumberTest.cpp
 *  \brief unit tests for VertexConditionNumberQualityMetric
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "VertexConditionNumberQualityMetric.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include "QualityMetricTester.hpp"

using namespace Mesquite;

class VertexConditionNumberTest : public CppUnit::TestFixture
{
private:  

  CPPUNIT_TEST_SUITE(VertexConditionNumberTest);

  CPPUNIT_TEST (test_supported_types);
  CPPUNIT_TEST (test_ideal_element_eval);
  CPPUNIT_TEST (test_ideal_element_grad);
  CPPUNIT_TEST (test_measures_quality);
  CPPUNIT_TEST (test_gradient_reflects_quality);
  CPPUNIT_TEST (test_inverted_elements);
  CPPUNIT_TEST (test_degenerate_elements);
  CPPUNIT_TEST (test_get_evaluations);
  CPPUNIT_TEST (test_get_element_indices);
  CPPUNIT_TEST (test_get_fixed_indices);
  CPPUNIT_TEST (test_eval_with_indices);
  CPPUNIT_TEST (test_location_invariant);
  CPPUNIT_TEST (test_orient_invariant);
  
  CPPUNIT_TEST_SUITE_END();
  
  VertexConditionNumberQualityMetric mMetric;
  QualityMetricTester tester;
  
public:

  VertexConditionNumberTest() 
    : tester(QualityMetricTester::NON_MIXED_FE)
    { }

  void test_supported_types()
    { tester.test_supported_element_types( &mMetric ); }
    
  void test_ideal_element_eval()
  {
    tester.test_evaluate_unit_tris_about_vertex( &mMetric, 1.0 );
    tester.test_evaluate_unit_quads_about_vertex( &mMetric, 1.0 );
    tester.test_evaluate_unit_hexes_about_vertex( &mMetric, 1.0 );

    tester.test_evaluate_unit_edge_element( &mMetric, TRIANGLE, 1.0 );
    tester.test_evaluate_unit_edge_element( &mMetric, QUADRILATERAL, 1.0 );
    tester.test_evaluate_unit_edge_element( &mMetric, TETRAHEDRON, 1.0 );
    tester.test_evaluate_unit_edge_element( &mMetric, HEXAHEDRON, 1.0 );
  }
  
  void test_ideal_element_grad()
    { tester.test_ideal_element_zero_vertex_gradient( &mMetric, false ); }
  
  void test_measures_quality()
    { tester.test_measures_vertex_quality( &mMetric ); }
  
  void test_gradient_reflects_quality()
    { tester.test_vertex_gradient_reflects_quality( &mMetric ); }
  
  void test_inverted_elements()
    { tester.test_evaluate_inverted_element( &mMetric, false ); }
    
  void test_degenerate_elements()
    { tester.test_evaluate_degenerate_element( &mMetric, false ); }
    
  void test_get_evaluations()
    { tester.test_get_vertex_evaluations( &mMetric ); }
    
  void test_get_element_indices()
    { tester.test_get_vertex_indices( &mMetric ); }
  
  void test_get_fixed_indices()
    { tester.test_get_indices_fixed( &mMetric ); }
  
  void test_eval_with_indices()
    { tester.compare_eval_and_eval_with_indices( &mMetric ); }
  
  void test_location_invariant()
  {
    tester.test_location_invariant( &mMetric );
    tester.test_grad_location_invariant( &mMetric );
  }
  
  void test_orient_invariant()
  {
    tester.test_orient_invariant( &mMetric );
    tester.test_grad_orient_invariant( &mMetric );
  }
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(VertexConditionNumberTest, "VertexConditionNumberTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(VertexConditionNumberTest, "Unit");
