/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "SundanceDoublingStepController.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceExprFieldWrapper.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"

namespace Sundance
{


bool DoublingStepController::run() const
{
  Tabs tab0(0);
  Tabs tab1;
  int verb = stepControl_.verbosity_;
  if (MPIComm::world().getRank()!=0) verb = 0;
  double t = stepControl_.tStart_;
  double p = stepControl_.stepOrder_;
  double tau = stepControl_.tau_;
  double tStop = stepControl_.tStop_;
  int maxSteps = stepControl_.maxSteps_;

  double safety = stepControl_.stepsizeReductionSafetyFactor_;


  double dt = stepControl_.initialStepsize_;

  double minStepsize = dt * stepControl_.minStepsizeFactor_;
  double maxStepsize = dt * stepControl_.maxStepsizeFactor_;
  bool atMinStepsize = false;
  bool atMaxStepsize = false;

  double minStepUsed = fabs(dt);
  double maxStepUsed = fabs(dt);

  PLAYA_MSG1(verb, tab0 << "=================================================="
    << endl
    << tab0 << "   starting time integration  " << endl
    << tab0 << "==================================================");
  PLAYA_MSG2(verb, tab1 << "Initial time: " << t);
  PLAYA_MSG2(verb, tab1 << "Final time: " << tStop);
  PLAYA_MSG2(verb, tab1 << "Initial stepsize: " << dt);
  PLAYA_MSG2(verb, tab1 << "Min allowed stepsize: " << minStepsize);
  PLAYA_MSG2(verb, tab1 << "Max allowed stepsize: " << maxStepsize);
  PLAYA_MSG2(verb, tab1 << "Step tolerance: " << tau);
  PLAYA_MSG2(verb, tab1 << "Method order: " << p);
  PLAYA_MSG2(verb, tab1 << "Max steps: " << maxSteps);

  Expr uCur = copyDiscreteFunction(prob_.uCur());
  Expr uTmp = copyDiscreteFunction(uCur);
  Expr uSolnOneFullStep = copyDiscreteFunction(uCur);
  Expr uSolnHalfStep = copyDiscreteFunction(uCur);
  Expr uSolnTwoHalfSteps = copyDiscreteFunction(uCur);

  int writeIndex = 1;
  double tNextWrite = outputControl_.writeInterval_;
  int writeVerb = outputControl_.verbosity_;
  double writeInterval = outputControl_.writeInterval_;

  /* Write the initial state */
  write(0, 0.0, uCur);
  

  int step = 0;
  int numShrink = 0;
  int numGrow = 0;
  int numSolves = 0;

  bool gotEvent = false;
  bool terminateOnDetection = false;
  if (eventHandler_.get()) terminateOnDetection 
    = eventHandler_->terminateOnDetection();
  PLAYA_MSG2(verb, tab1 << "Terminate on event detection: " 
    << terminateOnDetection);

  while (t < tStop && step < maxSteps)
  {
    Tabs tab2;
    PLAYA_MSG3(verb, tab2 << " step " << step << " time=" << t 
      << " to " << t + dt << " stepsize=" << dt);
    Tabs tab3;

    if (t + dt > tStop)
    {
      dt = tStop - t;
      PLAYA_MSG3(verb, tab3 << "timestep exceeds interval: reducing to dt=" << dt);
      PLAYA_MSG3(verb, tab3 << " step " << step << " is now time=" << t 
        << " to " << t + dt << " stepsize=" << dt);
    }

    PLAYA_MSG5(verb, tab3 << "doing full step");
    bool stepOK = prob_.step(t, uCur, t+dt, uSolnOneFullStep, solver_);
    PLAYA_MSG5(verb, tab3 << "doing first half step");
    stepOK = prob_.step(t, uCur, t+0.5*dt, uSolnHalfStep, solver_);
    PLAYA_MSG5(verb, tab3 << "doing second half step");
    stepOK = prob_.step(t+0.5*dt, uSolnHalfStep, t+dt, uSolnTwoHalfSteps, solver_);
    numSolves += 3;

    double err = compare_->diff(uSolnTwoHalfSteps, uSolnOneFullStep);

    PLAYA_MSG3(verb, tab3 << "error=" << err << " tau=" << tau);
    
    double dtNew = dt;
    if (err < tau) dtNew = dt*pow(tau/err, 1.0/(p+1.0));
    else dtNew = dt*pow(tau/err, 1.0/p);

    if (dtNew < dt && !atMinStepsize)
    {
      if (safety*dtNew <= minStepsize)
      {
        dt = minStepsize;
        atMinStepsize = true;
        PLAYA_MSG1(verb, tab2 << "WARNING: minimum stepsize reached");
      }
      else
      {
        dt = safety*dtNew;
      }
      atMaxStepsize = false;
      PLAYA_MSG3(verb, tab3 << "reducing step to dt=" << dt);
      numShrink++;
      continue;
    }

    if (stepHook_.get()) stepHook_->call(t+dt, uSolnTwoHalfSteps);

    if (eventHandler_.get())
    {
      gotEvent = eventHandler_->checkForEvent(t, uCur, t+0.5*dt, uSolnHalfStep);
      if (gotEvent && terminateOnDetection) break;

      gotEvent = eventHandler_->checkForEvent(t+0.5*dt, uSolnHalfStep,
        t+dt, uSolnTwoHalfSteps);
    }

    while ( tNextWrite >= t && tNextWrite < t+dt)
    {
      double tNode0 = t;
      double tNode1 = t+dt/2.0;
      double tNode2 = t+dt;
      
      double phi0 = (tNextWrite-tNode1)*(tNextWrite-tNode2)/(tNode0-tNode1)/(tNode0-tNode2);
      double phi1 = (tNextWrite-tNode0)*(tNextWrite-tNode2)/(tNode1-tNode0)/(tNode1-tNode2);
      double phi2 = (tNextWrite-tNode0)*(tNextWrite-tNode1)/(tNode2-tNode0)/(tNode2-tNode1);

      Vector<double> y0 = getDiscreteFunctionVector(uCur);
      Vector<double> y1 = getDiscreteFunctionVector(uSolnHalfStep);
      Vector<double> y2 = getDiscreteFunctionVector(uSolnTwoHalfSteps);

      Vector<double> yTmp = phi0*y0 + phi1*y1 + phi2*y2;
      setDiscreteFunctionVector(uTmp, yTmp);
      
      write(writeIndex, tNextWrite, uTmp);

      tNextWrite += writeInterval;
      writeIndex++;
    }

    minStepUsed = min(minStepUsed, fabs(dt));
    maxStepUsed = max(maxStepUsed, fabs(dt));

    updateDiscreteFunction(uSolnTwoHalfSteps, uCur);
    t += dt;
    step++;

    if (gotEvent && terminateOnDetection) break;

    if (t < tStop && dtNew > dt && !atMaxStepsize) /* increase timestep */
    {
      if (dtNew >= maxStepsize)
      {
        dt = maxStepsize;
        atMaxStepsize = true;
        PLAYA_MSG2(verb, tab2 << "WARNING: maximum stepsize reached");
      }
      else
      {
        dt = dtNew;
      }
      atMinStepsize = false;
      PLAYA_MSG3(verb, tab3 << "increasing step to dt=" << dt);
      numGrow++;
    }
  }

  

  PLAYA_MSG1(verb, tab0 << "=================================================================="
    << endl
    << tab0 << "   done time integration  ");
  PLAYA_MSG2(verb, tab0 << "------------------------------------------------------------------");
  PLAYA_MSG2(verb, tab1 << "Final time: " << t);
  PLAYA_MSG2(verb, tab1 << "Num steps: " << step);
  PLAYA_MSG2(verb, tab1 << "Min stepsize used: " << minStepUsed
    << " (would need " << fabs(t-stepControl_.tStart_)/minStepUsed 
    << " steps of this size)");
  PLAYA_MSG2(verb, tab1 << "Max stepsize used: " << maxStepUsed);
  PLAYA_MSG2(verb, tab1 << "Num nonlinear solves: " << numSolves);
  PLAYA_MSG2(verb, tab1 << "Num step reductions: " << numShrink);
  PLAYA_MSG2(verb, tab1 << "Num step increases: " << numGrow);
  if (terminateOnDetection && gotEvent)
  {
    PLAYA_MSG2(verb, tab1 << "Terminated by event detection at time="
      << eventHandler_->eventTime());
  }
  PLAYA_MSG1(verb, tab0 << "==================================================================");

}


void DoublingStepController::write(int index, double t, const Expr& u) const
{
  Tabs tab(0);
  string name = outputControl_.filename_ + "-" + Teuchos::toString(index);

  FieldWriter w = outputControl_.wf_.createWriter(name);

  PLAYA_MSG1(outputControl_.verbosity_, tab << "writing to " << name);

  Mesh mesh = getDiscreteFunctionMesh(u);
  w.addMesh(mesh);
  for (int i=0; i<u.size(); i++)
  {
    w.addField("output[" + Teuchos::toString(i) + "]", 
      new ExprFieldWrapper(u[i]));
  }
  w.write();
  
}


}
