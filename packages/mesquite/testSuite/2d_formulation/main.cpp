/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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

    (2007) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file main.cpp
 *  \brief Try examples from "Formulation of a Target-Matrix Paradigm
 *         for Mesh Optimization", Patrick Knupp.
 *  \author Jason Kraftcheck 
 */
 
#define USE_GLOBAL_PATCH

#include "Mesquite.hpp"

#include "PMeanPTemplate.hpp"
#include "AffineMapMetric.hpp"
#include "ConjugateGradient.hpp"
#include "TerminationCriterion.hpp"
#include "ElementPMeanP.hpp"
#include "MsqError.hpp"
#include "LinearFunctionSet.hpp"
#include "UnitWeight.hpp"
#include "TSquared2D.hpp"
#include "SamplePoints.hpp"
#include "MeshImpl.hpp"
#include "PlanarDomain.hpp"
#include "InstructionQueue.hpp"
#include "IdealTargetCalculator.hpp"
#include "IdentityTarget.hpp"
#include "MetricWeight.hpp"
#include "InverseMetricWeight.hpp"
#include "TargetWriter.hpp"
#include "WeightReader.hpp"

#include <iostream>
#include <stdlib.h>

using namespace Mesquite;
using namespace std;

const double epsilon = 2e-2;
const bool write_results = true;

#define CHKERR(A) if (A) { cerr << (A) << endl; exit(1); }

enum Grouping { SAMPLE, ELEMENT, QUADRANT, HALF };
enum Weight { UNIT, METRIC, INV_METRIC };

void run_test( Grouping grouping, int of_power, Weight w, const string filename )
{
  MsqError err;
  
  UnitWeight wc;
  IdentityTarget target;
  TSquared2D target_metric;
  AffineMapMetric qual_metric( &target, &wc, &target_metric, NULL );
  ElementPMeanP elem_metric( of_power, &qual_metric );
  QualityMetric* qm_ptr = (grouping == ELEMENT) ? (QualityMetric*)&elem_metric : (QualityMetric*)&qual_metric;

  PMeanPTemplate OF( of_power, qm_ptr );
  ConjugateGradient solver( &OF, true );
  TerminationCriterion tc;
  TerminationCriterion itc;
  tc.add_criterion_type_with_double( TerminationCriterion::VERTEX_MOVEMENT_ABSOLUTE, 1e-4, err );
  itc.add_criterion_type_with_int( TerminationCriterion::NUMBER_OF_ITERATES, 2, err );
#ifdef USE_GLOBAL_PATCH
  solver.use_global_patch();
  solver.set_inner_termination_criterion( &tc );
#else
  solver.use_element_on_vertex_patch();
  solver.set_inner_termination_criterion( &itc );
  solver.set_outer_termination_criterion( &tc );
#endif
  
  MeshImpl mesh, expected_mesh;
  mesh.read_vtk( SRCDIR "/initial.vtk", err ); CHKERR(err)
//  expected_mesh.read_vtk( (filename + ".vtk").c_str(), err ); CHKERR(err)
  
  PlanarDomain plane( PlanarDomain::XY );
  LinearFunctionSet maps;

  MetricWeight mw( &qual_metric );
  InverseMetricWeight imw( &qual_metric );
  WeightReader reader;
  if (w == METRIC) {
    TargetWriter writer( qual_metric.get_sample_points(), 0, &mw );
    InstructionQueue tq;
    tq.add_target_calculator( &writer, err );
    tq.run_instructions( &mesh, &plane, &maps, err ); CHKERR(err);
    qual_metric.set_weight_calculator( &reader );
  }
  else if (w == INV_METRIC) {
    TargetWriter writer( qual_metric.get_sample_points(), 0, &imw );
    InstructionQueue tq;
    tq.add_target_calculator( &writer, err );
    tq.run_instructions( &mesh, &plane, &maps, err ); CHKERR(err);
    qual_metric.set_weight_calculator( &reader );
  }

  InstructionQueue q;
  q.set_master_quality_improver( &solver, err );
  q.run_instructions( &mesh, &plane, &maps, err ); CHKERR(err)
/*  
  vector<Mesh::VertexHandle> vemain.cpprts;
  vector<MsqVertex> mesh_coords, expected_coords;
  mesh.get_all_vertices( verts, err ); CHKERR(err)
  mesh_coords.resize(verts.size());
  mesh.vertices_get_coordinates( &verts[0], &mesh_coords[0], verts.size(), err ); CHKERR(err)
  expected_mesh.get_all_vertices( verts, err ); CHKERR(err)
  expected_coords.resize(verts.size());
  expected_mesh.vertices_get_coordinates( &verts[0], &expected_coords[0], verts.size(), err ); CHKERR(err)
  if (expected_coords.size() != mesh_coords.size()) {
    cerr << "Invlid expected mesh.  Vertex count doesn't match" << endl;
    exit(1);
  }
  
  unsigned error_count = 0;
  for (size_t i = 0; i < mesh_coords.size(); ++i) 
    if ((expected_coords[i] - mesh_coords[i]).length_squared() > epsilon*epsilon)
      ++error_count;

  if (!error_count) 
    cout << filename << " : SUCCESS" << endl;
  else
    cout << filename << " : FAILURE (" << error_count 
         << " vertices differ by more than " << epsilon << ")" << endl;
*/
  if (write_results)
    mesh.write_vtk( (filename + ".results.vtk").c_str(), err ); CHKERR(err)
}

int main()
{
  run_test( SAMPLE, 1, UNIT, "1-1" );
  run_test( SAMPLE, 2, UNIT, "1-2" );
  run_test( SAMPLE, 4, UNIT, "1-4" );
  run_test( SAMPLE, 8, UNIT, "1-8" );
  
  run_test(  SAMPLE, 1, UNIT, "2-NW" );
  run_test( ELEMENT, 1, UNIT, "2-NE" );
  
  run_test( SAMPLE, 1,       UNIT, "3-Left"  );
  run_test( SAMPLE, 1,     METRIC, "3-Mid"   );
  run_test( SAMPLE, 1, INV_METRIC, "3-Right" );
}
