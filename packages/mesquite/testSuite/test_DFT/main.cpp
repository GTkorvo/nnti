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
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 13-Apr-04 at 10:57:52
//  LAST-MOD: 21-May-04 at 13:45:21 by Thomas Leurent
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


#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MeshMunson.hpp"
#include "MesquiteError.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorythms
#include "I_DFT.hpp"
#include "sI_DFT.hpp"
#include "RI_DFT.hpp"
#include "sRI_DFT.hpp"
#include "ConcreteTargetCalculators.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
#include "MsqMessage.hpp"
using namespace Mesquite;

using std::cout;
using std::endl;

#undef __FUNC__
#define __FUNC__ "main"
int main(int argc, char* argv[])
{
  Mesquite::MsqError err;
  char file_name[256];
  double OF_value = 1.;
  
  // command line arguments
  if (argc==1 || argc>3)
    cout << "meshfile name needed as argument.\n"
      "objective function value optional as 2nd argument.\n" << endl;
  else if (argc==2) {
    cout << " given 1 command line argument.\n";
    strcpy(file_name, argv[1]);
  } else if (argc==3) {
    cout << " given 2 command line arguments.\n";
    strcpy(file_name, argv[1]);
    OF_value = atof(argv[2]);
  }
  
  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  mesh->read_vtk(file_name, err);
  
  // initialises a MeshSet object
  MeshSet mesh_set1;
  mesh_set1.add_mesh(mesh, err); MSQ_CHKERR(err);

  // creates an intruction queue
  InstructionQueue queue1;

  // creates a mean ratio quality metric ...
  I_DFT mean_ratio;
//  mean_ratio.set_gradient_type(QualityMetric::NUMERICAL_GRADIENT);
//   mean_ratio.set_hessian_type(QualityMetric::NUMERICAL_HESSIAN);
   mean_ratio.set_averaging_method(QualityMetric::LINEAR, err); MSQ_CHKERR(err);

  // creates a target calculator
//  DefaultTargetCalculator target;

  Mesquite::MeshImpl *ref_mesh = new Mesquite::MeshImpl;
  ref_mesh->read_vtk("../../meshFiles/3D/VTK/cube_hex_rm.vtk", err);
  MeshSet ref_mesh_set;
  ref_mesh_set.add_mesh(ref_mesh, err); MSQ_CHKERR(err);
  ReferenceMeshTargetCalculator target(&ref_mesh_set);

  Mesquite::MeshImpl *ref_mesh2 = new Mesquite::MeshImpl;
  ref_mesh2->read_vtk("../../meshFiles/3D/VTK/cube_hex_rm.vtk", err);
  MeshSet ref_mesh2_set;
  ref_mesh2_set.add_mesh(ref_mesh2, err); MSQ_CHKERR(err);
  ReferenceMeshTargetCalculator assessor_target(&ref_mesh2_set);
   
  // ... and builds an objective function with it
  LPtoPTemplate* obj_func = new LPtoPTemplate(&mean_ratio, 1, err);
  obj_func->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
  
  // creates the steepest descentfeas newt optimization procedures
//  ConjugateGradient* pass1 = new ConjugateGradient( obj_func, err );
  FeasibleNewton* pass1 = new FeasibleNewton( obj_func );
  pass1->set_target_calculator(&target, err); MSQ_CHKERR(err);
  pass1->set_patch_type(PatchData::GLOBAL_PATCH, err);
  
  QualityAssessor stop_qa(&mean_ratio,QualityAssessor::AVERAGE);
  stop_qa.set_target_calculator(&assessor_target, err); MSQ_CHKERR(err);
  
  // **************Set stopping criterion****************
  TerminationCriterion tc_inner;
  tc_inner.add_criterion_type_with_double(
           TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE, OF_value, err);
  TerminationCriterion tc_outer;
  tc_outer.add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  pass1->set_inner_termination_criterion(&tc_inner);
  pass1->set_outer_termination_criterion(&tc_outer);

  // sets a culling method on the first QualityImprover
  pass1->add_culling_method(PatchData::NO_BOUNDARY_VTX);
  
  queue1.add_quality_assessor(&stop_qa, err); MSQ_CHKERR(err);
   
  // adds 1 pass of pass1 to mesh_set1
  queue1.set_master_quality_improver(pass1, err); MSQ_CHKERR(err);
  
  queue1.add_quality_assessor(&stop_qa, err); MSQ_CHKERR(err);

mesh->write_vtk("original_mesh",err); MSQ_CHKERR(err);
  
  // launches optimization on mesh_set1
 queue1.run_instructions(mesh_set1, err); MSQ_CHKERR(err);
  
mesh->write_vtk("smoothed_mesh", err); MSQ_CHKERR(err);
  PRINT_TIMING_DIAGNOSTICS();
}
