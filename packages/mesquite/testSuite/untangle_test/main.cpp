/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 23-Jul-03 at 18:08:13 by Thomas Leurent
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#include "meshfiles.h"

#include "MeshImpl.hpp"
#include "MsqTimer.hpp"
#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorithms
#include "Randomize.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"
#include "ConjugateGradient.hpp"
#include "PlanarDomain.hpp"

#include "UntangleWrapper.hpp"

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>

using namespace Mesquite;

#define VTK_2D_DIR MESH_FILES_DIR "2D/VTK/"

// This was the original 'main' code before tests of UntangleWrapper were added
int old_untangle_beta_test();

// Test untangle wrapper
// Assumes all meshes lie in a plane for which the normal is [0,0,1].
int uwt( bool skip,
         UntangleWrapper::UntangleMetric metric,
         const char* input_file,
         int expected_number_of_remaining_inverted_elems,
         bool flip_domain = false );

bool brief_output = false;
bool write_output = false;
double mu_sigma = -1;
double beta = -1;

void usage( const char* argv0 ) {
  std::cerr << "Usage: " << argv0 << " [-<flags>] [-c <sigma>] [-b <beta>]" << std::endl
            << "       " << argv0 << " -h" << std::endl;
}
void help( const char* argv0 ) {
  std::cout << "Usage: " << argv0 << " [-<flags>] [-c <sigma>] [-b <beta>]" << std::endl
            << "Flags: -q : brief output" << std::endl
            << "       -w : write result meshes" << std::endl
            << "       -O : skip legacy untangle beta test" << std::endl
            << "       -B : skip tests using untangle beta target metric" << std::endl
            << "       -Z : skip tests using size untangle target metric" << std::endl
            << "       -P : skip tests using shapesize untangle target metric" << std::endl
            << "       -H : skip tests using 'tangled_horse1.vtk' as input" << std::endl
            << "       -Q : skip tests using 'hole_in_square_tanglap.vtk' as input" << std::endl
            << "       -I : skip tests using 'inverted-hole-2.vtk' as input" << std::endl
            << "       -S : skip tests using 'shest_grid32.vtk' as input" << std::endl
            << "       -c : specify sigma value for untangle metrics" << std::endl
            << "       -b : specify beta value for untangle beta metric" << std::endl
            << std::endl;
}

int main( int argc, char* argv[] )
{
  bool skip_old = false;
  bool skip_beta = false;
  bool skip_size = false;
  bool skip_shape = false;
  bool skip_horse = false;
  bool skip_hole = false;
  bool skip_invrt = false;
  bool skip_shest = false;
  std::list<double*> expected;

  for (int i = 1; i < argc; ++i) {
    if (!expected.empty()) {
      char* endptr;
      *expected.front() = strtod( argv[i], &endptr );
      if (*endptr || *expected.front() <= 0) {
        std::cerr << "Expected positive number, found \"" << argv[i] << '"' << std::endl;
        usage(argv[0]);
        return 1;
      }
      expected.pop_front();
    }
    else if (argv[i][0] == '-' && argv[i][1]) {
      for (int j = 1; argv[i][j]; ++j) {
        switch (argv[i][j]) {
          case 'q': brief_output = true; break;
          case 'w': write_output = true; break;
          case 'O': skip_old = true; break;
          case 'B': skip_beta = true; break;
          case 'Z': skip_size = true; break;
          case 'P': skip_shape = true; break;
          case 'H': skip_horse = true; break;
          case 'Q': skip_hole = true; break;
          case 'I': skip_invrt = true; break;
          case 'S': skip_shest = true; break;
          case 'c': expected.push_back(&mu_sigma); break;
          case 'b': expected.push_back(&beta); break;
          case 'h': help(argv[0]); return 0;
          default:
            std::cerr << "Invalid flag: -" << argv[i][j] << std::endl;
            usage(argv[0]);
            return 1;
        }
      }
    }
    else {
      std::cerr << "Unexpected argument: \"" << argv[i] << '"' << std::endl;
      usage(argv[0]);
      return 1;
    }
  }

  int result = 0;
  if (!skip_old)
    result = old_untangle_beta_test();
  
  result += uwt( skip_beta||skip_horse, UntangleWrapper::BETA, "tangled_horse1.vtk",         0 );
  result += uwt( skip_beta||skip_hole , UntangleWrapper::BETA, "hole_in_square_tanglap.vtk", 0, true );
  result += uwt( skip_beta||skip_invrt, UntangleWrapper::BETA, "inverted-hole-2.vtk",        0 );
  result += uwt( skip_beta||skip_shest, UntangleWrapper::BETA, "shest_grid32.vtk",           0 );

  result += uwt( skip_size||skip_horse, UntangleWrapper::SIZE, "tangled_horse1.vtk",         0 );
  result += uwt( skip_size||skip_hole , UntangleWrapper::SIZE, "hole_in_square_tanglap.vtk", 6, true );
  result += uwt( skip_size||skip_invrt, UntangleWrapper::SIZE, "inverted-hole-2.vtk",        0  );
  result += uwt( skip_size||skip_shest, UntangleWrapper::SIZE, "shest_grid32.vtk",           0 );

  result += uwt( skip_shape||skip_horse, UntangleWrapper::SHAPESIZE, "tangled_horse1.vtk",         0 );
  result += uwt( skip_shape||skip_hole , UntangleWrapper::SHAPESIZE, "hole_in_square_tanglap.vtk", 0, true );
  result += uwt( skip_shape||skip_invrt, UntangleWrapper::SHAPESIZE, "inverted-hole-2.vtk",        8  );
  result += uwt( skip_shape||skip_shest, UntangleWrapper::SHAPESIZE, "shest_grid32.vtk",           0 );
  
  return result;
}

int old_untangle_beta_test()
{
  Mesquite::MeshImpl mesh;
  MsqPrintError err(cout);
  mesh.read_vtk(VTK_2D_DIR "tangled_quad.vtk", err);
  if (err) return 1;
  
  // Set Domain Constraint
  Vector3D pnt(0,0,0);
  Vector3D s_norm(0,0,1);
  PlanarDomain msq_geom(s_norm, pnt);
                                                                              
    // creates an intruction queue
  InstructionQueue queue1;
  
    // creates a mean ratio quality metric ...
  ConditionNumberQualityMetric shape_metric;
  UntangleBetaQualityMetric untangle(2);
  Randomize pass0(.05);
    // ... and builds an objective function with it
    //LInfTemplate* obj_func = new LInfTemplate(shape_metric);
  LInfTemplate obj_func(&untangle);
  LPtoPTemplate obj_func2(&shape_metric, 2, err);
  if (err) return 1;
    // creates the steepest descent optimization procedures
  ConjugateGradient pass1( &obj_func, err );
  if (err) return 1;
  
    //SteepestDescent* pass2 = new SteepestDescent( obj_func2 );
  ConjugateGradient pass2( &obj_func2, err );
  if (err) return 1;
  pass2.use_element_on_vertex_patch();
  if (err) return 1;
  pass2.use_global_patch();
  if (err) return 1;
  QualityAssessor stop_qa=QualityAssessor(&shape_metric);
  QualityAssessor stop_qa2=QualityAssessor(&shape_metric);
  if (brief_output) {
    stop_qa.disable_printing_results();
    stop_qa2.disable_printing_results();
  }
  
  stop_qa.add_quality_assessment(&untangle);
    // **************Set stopping criterion**************
    //untangle beta should be 0 when untangled
  TerminationCriterion sc1;
  sc1.add_relative_quality_improvement( 0.000001 );
  TerminationCriterion sc3;
  sc3.add_iteration_limit( 10 );
  TerminationCriterion sc_rand;
  sc_rand.add_iteration_limit( 1 );
  
    //StoppingCriterion sc1(&stop_qa,-1.0,.0000001);
    //StoppingCriterion sc3(&stop_qa2,.9,1.00000001);
    //StoppingCriterion sc2(StoppingCriterion::NUMBER_OF_PASSES,10);
    //StoppingCriterion sc_rand(StoppingCriterion::NUMBER_OF_PASSES,1);
    //either until untangled or 10 iterations
  pass0.set_outer_termination_criterion(&sc_rand);
  pass1.set_outer_termination_criterion(&sc1);
  pass2.set_inner_termination_criterion(&sc3);
  
    // adds 1 pass of pass1 to mesh_set1
  queue1.add_quality_assessor(&stop_qa,err); 
  if (err) return 1;
    //queue1.add_preconditioner(pass0,err);MSQ_CHKERR(err);
    //queue1.add_preconditioner(pass1,err);MSQ_CHKERR(err);
    //queue1.set_master_quality_improver(pass2, err); MSQ_CHKERR(err);
  queue1.set_master_quality_improver(&pass1, err);
  if (err) return 1;
  queue1.add_quality_assessor(&stop_qa2,err);
  if (err) return 1;
  if (write_output)
    mesh.write_vtk("original_mesh.vtk", err);
  if (err) return 1;
  
    // launches optimization on mesh_set1
  queue1.run_instructions(&mesh, &msq_geom, err);
  if (err) return 1;
  
  if (write_output)
    mesh.write_vtk("smoothed_mesh.vtk", err); 
  if (err) return 1;
  
  if (!brief_output)
    print_timing_diagnostics(cout);
  return 0;
}

const char* tostr( UntangleWrapper::UntangleMetric m )
{
  static const char BETA[] = "BETA";
  static const char SIZE[] = "SIZE";
  static const char SHAPESIZE[] = "SHAPESIZE";
  switch (m) {
    case UntangleWrapper::BETA:      return BETA;
    case UntangleWrapper::SIZE:      return SIZE;
    case UntangleWrapper::SHAPESIZE: return SHAPESIZE;
  }
  return 0;
}

int uwt( bool skip,
         UntangleWrapper::UntangleMetric metric,
         const char* input_file_base,
         int expected,
         bool flip_domain )
{
  if (skip)
    return 0;

  if (!brief_output)
    std::cout << std::endl << "**********************************************" << std::endl;
  std::cout << "Running \"" << input_file_base << "\" for " << tostr(metric) << std::endl;
  if (!brief_output)
    std::cout << "**********************************************" << std::endl << std::endl;

    // get mesh
  MsqError err;
  MeshImpl mesh;
  std::string input_file( VTK_2D_DIR );
  input_file += input_file_base;
  mesh.read_vtk( input_file.c_str(), err );
  if (err) {
    std::cerr << err << std::endl;
    std::cerr << "ERROR: " << input_file << " : failed to read file" << std::endl;
    return 1;
  }
    // get domain
  std::vector<Mesh::VertexHandle> verts;
  mesh.get_all_vertices( verts, err );
  if (err || verts.empty()) abort();
  MsqVertex coords;
  mesh.vertices_get_coordinates( arrptr(verts), &coords, 1, err );
  if (err) abort();
  Vector3D norm(0,0,flip_domain ? -1 : 1);
  PlanarDomain domain( norm, coords );
    // run wrapper
  UntangleWrapper wrapper( metric );
  wrapper.set_vertex_movement_limit_factor( 0.01 );
  double constant = (metric == UntangleWrapper::BETA) ? beta : mu_sigma;
  if (constant > 0)
    wrapper.set_metric_constant( constant );
  if (brief_output)
    wrapper.quality_assessor().disable_printing_results();
  wrapper.run_instructions( &mesh, &domain, err );
  if (err) {
    std::cerr << err << std::endl;
    std::cerr << "ERROR: optimization failed" << std::endl;
    return 1;
  }
    // write output file
  if (write_output) {
    std::string result_file(tostr(metric));
    result_file += "-";
    result_file += input_file_base;
    mesh.write_vtk( result_file.c_str(), err );
    if (err) {
      std::cerr << err << std::endl;
      std::cerr << "ERROR: " << result_file << " : failed to write file" << std::endl;
      err.clear();
    }
    else {
      std::cerr << "Wrote file: " << result_file << std::endl;
    }
  }
  
    // test number of inverted elements
  int count, junk;
  wrapper.quality_assessor().get_inverted_element_count( count, junk, err );
  if (err) abort();
  if (count < expected) {
    std::cout << "WARNING: expected " << expected 
              << " inverted elements but finished with only " 
              << count << std::endl
              << "Test needs to be updated?" << std::endl << std::endl;
    return 0;
  }
  else if (count == expected) {
    std::cout << "Completed with " << count << " inverted elements remaining" 
              << std::endl << std::endl;
    return 0;
  }
  else {
    std::cerr << "ERROR: expected " << expected 
              << " inverted elements but finished with " 
              << count << std::endl << std::endl;
    return 1;
  }
}
