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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file ShapeImprovementWrapper.cpp

Member functions of the Mesquite::ShapeImprovementWrapper class

  \author Michael Brewer
  \date   June 6, 2003
 */

#include "InstructionQueue.hpp"
#include "ShapeImprovementWrapper.hpp"
#include "MeshSet.hpp"
#include "MsqTimer.hpp"

using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "ShapeImprovementWrapper::ShapeImprovementWrapper"
/*! The consturctor allows for two values.  The first is a 
  time bound (in seconds) used as a termination criterion.  If
  this value is non-positive, no time bound will be set.
  By default, the value is set to zero and no time bound
  is used.  The second value is the tolerance for the gradient
  norm termination criteria.  The default value is 1.e-6.*/
ShapeImprovementWrapper::ShapeImprovementWrapper(double cpu_time,
                                                 double grad_norm) {

    //arbitrarily chosen variables
  untBeta=1.e-8;
  successiveEps=1.e-4;
  
  
  
  if(cpu_time>0.0){
    timerNeeded=true;
  }
  else{
    timerNeeded=false;
  }
  maxTime=cpu_time;
  
  MsqError err;
  untangleMetric = new UntangleBetaQualityMetric(untBeta);
  untangleFunc =  new LPtoPTemplate(untangleMetric, 2, err);
  untangleFunc->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
  untangleGlobal = new ConjugateGradient(untangleFunc,err);
  untangleGlobal->set_patch_type(PatchData::GLOBAL_PATCH, err,1 ,1);
  
  untangleGlobalInner = new TerminationCriterion();
  untangleGlobalOuter = new TerminationCriterion();
  
  untangleGlobalInner->add_criterion_type_with_double(TerminationCriterion::QUALITY_IMPROVEMENT_ABSOLUTE,0.0,err);
  untangleGlobalInner->add_criterion_type_with_double(TerminationCriterion::SUCCESSIVE_IMPROVEMENTS_ABSOLUTE,successiveEps,err);
  untangleGlobalOuter->add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
  
  meanRatio = new MeanRatioQualityMetric;
  meanRatio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
  meanRatio->set_hessian_type(QualityMetric::ANALYTICAL_HESSIAN);
  meanRatio->set_averaging_method(QualityMetric::LINEAR,err);
    // creates the l_2 squared objective function
  objFunc = new LPtoPTemplate(meanRatio, 2, err);
  objFunc->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
    //creates a FeasibleNewtone improver
  feasNewt = new FeasibleNewton(objFunc);
  feasNewt->set_patch_type(PatchData::GLOBAL_PATCH, err,1 ,1);
  mQA = new QualityAssessor(meanRatio,QualityAssessor::MAXIMUM);
  mQA->add_quality_assessment(meanRatio, QualityAssessor::MINIMUM,err);
  mQA->add_quality_assessment(meanRatio, QualityAssessor::AVERAGE,err);
  mQA->add_quality_assessment(meanRatio, QualityAssessor::RMS,err);   
        //**************Set stopping criterion*e***************
  termInner = new TerminationCriterion();
  termOuter = new TerminationCriterion();
  termInner->add_criterion_type_with_double(TerminationCriterion::GRADIENT_L2_NORM_ABSOLUTE,grad_norm,err);
  termInner->add_criterion_type_with_double(TerminationCriterion::SUCCESSIVE_IMPROVEMENTS_RELATIVE,successiveEps,err);
  termOuter->add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err);
    // sets a culling method on the first QualityImprover
  untangleGlobal->add_culling_method(PatchData::NO_BOUNDARY_VTX);
  untangleGlobal->set_inner_termination_criterion(untangleGlobalInner);
  untangleGlobal->set_outer_termination_criterion(untangleGlobalOuter);
    // sets a culling method on the second QualityImprover
    //untangleLocal->add_culling_method(PatchData::NO_BOUNDARY_VTX);
    //untangleLocal->set_inner_termination_criterion(untangleLocalInner);
    //untangleLocal->set_outer_termination_criterion(untangleLocalOuter);
    // sets a culling method on the third QualityImprover
  feasNewt->add_culling_method(PatchData::NO_BOUNDARY_VTX);
  feasNewt->set_inner_termination_criterion(termInner);
  feasNewt->set_outer_termination_criterion(termOuter);
      
}


#undef __FUNC__
#define __FUNC__ "ShapeImprovementWrapper::~ShapeImprovementWrapper"
ShapeImprovementWrapper::~ShapeImprovementWrapper()
{
  delete untangleMetric;
  delete untangleFunc;
  delete untangleGlobal;
  delete untangleGlobalInner;
  delete untangleGlobalOuter;
      
  delete meanRatio;
  delete objFunc;
  delete feasNewt;
  delete mQA;
  delete termInner;
  delete termOuter;
}


#undef __FUNC__
#define __FUNC__ "ShapeImprovementWrapper::run_instructions"
/*!Run instructions first calls the global untangler.  If the
  resulting mesh is tangled after that pre-conditioning step,
  The mesh is iteratively smoothed with a local and then global
  untangler until the mesh is untangled or until a certain time
  constraint has been exceeded.  If the mesh was successfully
  untangled and there is still time remaining, a mean ratio
  shape improvement is then performed.*/
void ShapeImprovementWrapper::run_instructions(MeshSet &ms, MsqError &err)
{
    //a timer to keep track of the amount of time spent in this wrapper
  Timer totalTimer;
    //time remaining keeps a track of how much time is left before the
    //wrapper must terminate.  If the wrapper is set to terminate on
    //a time constraint, time_remaining will always be 1.0
  double time_remaining=1.0;
  mQA->loop_over_mesh(ms, err);
    //if using a time constraint set the termination criteria.
  if(timerNeeded){
    time_remaining=maxTime;
    untangleGlobalInner->add_criterion_type_with_double(TerminationCriterion::CPU_TIME,time_remaining,err);
  }
    //global untangler
  untangleGlobal->loop_over_mesh(ms, err);
  if(timerNeeded)
    time_remaining=maxTime-totalTimer.since_birth();
  double func_val=untangleGlobalInner->get_current_function_value();

  mQA->loop_over_mesh(ms, err);
  if(timerNeeded)
    time_remaining=maxTime-totalTimer.since_birth();
    //if all the time constraint has been exceeded, notify the user that
    //the shape improvement has not been performed.
  if(time_remaining<=0){
    Message::print_info("\nOptimization is terminating without perfoming shape improvement");
    Message::print_info("\n Untangle Function Value is %f",func_val);
  }
    //otherwise, perform the shape improvement.
  else{
    if(timerNeeded)
      termInner->add_criterion_type_with_double(TerminationCriterion::CPU_TIME,time_remaining,err);
    feasNewt->loop_over_mesh(ms, err);
    mQA->loop_over_mesh(ms, err);
  } 
}

  
