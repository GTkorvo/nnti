/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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


/** \file BCDTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "LPtoPTemplate.hpp"
#include "PowerMeanP.hpp"
#include "IdealWeightInverseMeanRatio.hpp"
#include "InstructionQueue.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
#include "MeshImpl.hpp"
#include "QualityAssessor.hpp"
#include "TerminationCriterion.hpp"

#include "UnitUtil.hpp"
#include "meshfiles.h"
#include "cppunit/extensions/HelperMacros.h"
#include <iostream>
#include <vector>
#include <string>

using namespace Mesquite;
using namespace std;

const char HEX_MESH[] = MESH_FILES_DIR "3D/VTK/1000hex-block-internal-bias.vtk";
const char TET_MESH[] = MESH_FILES_DIR "3D/VTK/tire.vtk";
typedef FeasibleNewton SolverType;

class BCDTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(BCDTest);

  CPPUNIT_TEST (test_lp_to_p_hex);
  CPPUNIT_TEST (test_p_mean_p_hex);
  CPPUNIT_TEST (test_lp_to_p_tet);
  CPPUNIT_TEST (test_p_mean_p_tet);

  CPPUNIT_TEST_SUITE_END();

  IdealWeightInverseMeanRatio mMetric;

  void compare_bcd( ObjectiveFunction* of, string name, const char* file );
  
public:

  void test_lp_to_p_hex() {
    LPtoPTemplate OF( 1, &mMetric );
    compare_bcd( &OF, "LPtoP-hex", HEX_MESH );
  }
  
  void test_p_mean_p_hex() {
    PowerMeanP OF( 1.0, &mMetric );
    compare_bcd( &OF, "PMeanP-hex", HEX_MESH );
  }

  void test_lp_to_p_tet() {
    LPtoPTemplate OF( 1, &mMetric );
    compare_bcd( &OF, "LPtoP-tet", TET_MESH );
  }
  
  void test_p_mean_p_tet() {
    PowerMeanP OF( 1.0, &mMetric );
    compare_bcd( &OF, "PMeanP-tet", TET_MESH );
  }
  
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(BCDTest, "BCDTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(BCDTest, "Regression");

void BCDTest::compare_bcd( ObjectiveFunction* OF, string name, const char* mesh_file )
{
  MsqPrintError err(cout);
  size_t i;
  vector<MsqVertex> initial_coords, global_coords, bcd_coords;
  vector<Mesh::VertexHandle> vertex_list;
  
    // set up a smoother
  TerminationCriterion iterations, vertex_movement;
  iterations.add_criterion_type_with_int( TerminationCriterion::NUMBER_OF_ITERATES, 2, err );
  vertex_movement.add_criterion_type_with_double( TerminationCriterion::VERTEX_MOVEMENT_ABSOLUTE, 1e-3, err );

  SolverType global_solver( OF, true );
  SolverType bcd_solver( OF, false );
  global_solver.use_global_patch();
  bcd_solver.use_element_on_vertex_patch();
  global_solver.set_inner_termination_criterion( &vertex_movement );
  bcd_solver.set_inner_termination_criterion( &iterations );
  bcd_solver.set_outer_termination_criterion( &vertex_movement );

  QualityAssessor qa;
  const int NO_HIST = QualityAssessor::ALL_MEASURES & ~QualityAssessor::HISTOGRAM;
  qa.add_quality_assessment( &mMetric, NO_HIST, err );

  InstructionQueue global_q, bcd_q;
  global_q.add_quality_assessor( &qa, err );
  global_q.set_master_quality_improver( &global_solver, err );
  global_q.add_quality_assessor( &qa, err );
  bcd_q.set_master_quality_improver( &bcd_solver, err );
  bcd_q.add_quality_assessor( &qa, err );
  
    // read mesh
  MeshImpl mesh;
  mesh.read_vtk( mesh_file, err ); ASSERT_NO_ERROR(err);
  mesh.get_all_vertices( vertex_list, err ); ASSERT_NO_ERROR(err);
  initial_coords.resize( vertex_list.size() );
  mesh.vertices_get_coordinates( &vertex_list[0], &initial_coords[0], vertex_list.size(), err );
  ASSERT_NO_ERROR(err);
  
    // run global smoother
  global_q.run_instructions( &mesh, err ); 
  ASSERT_NO_ERROR(err);
  mesh.write_vtk( (name + "-gbl.vtk").c_str(), err );
  global_coords.resize( vertex_list.size() );
  mesh.vertices_get_coordinates( &vertex_list[0], &global_coords[0], vertex_list.size(), err );
  ASSERT_NO_ERROR(err);
  
    // restore initial vertex positions
  for (i = 0; i < vertex_list.size(); ++i) {
    mesh.vertex_set_coordinates( vertex_list[i], initial_coords[i], err );
    ASSERT_NO_ERROR(err);
  }
  
    // run local smoother
  bcd_q.run_instructions( &mesh, err );
  ASSERT_NO_ERROR(err);
  mesh.write_vtk( (name + "-bcd.vtk").c_str(), err );
  bcd_coords.resize( vertex_list.size() );
  mesh.vertices_get_coordinates( &vertex_list[0], &bcd_coords[0], vertex_list.size(), err );
  ASSERT_NO_ERROR(err);
  
    // compare results
  for (i = 0; i < bcd_coords.size(); ++i)
    CPPUNIT_ASSERT_VECTORS_EQUAL( global_coords[i], bcd_coords[i], 1e-2 );
}
