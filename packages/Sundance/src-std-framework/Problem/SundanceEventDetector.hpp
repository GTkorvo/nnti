/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef SUNDANCE_EVENT_DETECTOR_H
#define SUNDANCE_EVENT_DETECTOR_H


#include "SundanceDefs.hpp"

namespace Sundance
{
class Expr;

/** */
class EventDetectorBase
{
public:
  /** */
  EventDetectorBase() {}
  
  /** */
  virtual bool terminateOnDetection() const {return false;}

  /** */
  virtual bool checkForEvent(
    const double& t1, const Expr& u1,
    const double& t2, const Expr& u2) = 0 ;

  /** */
  virtual double eventTime() const = 0 ;

  /** */
  virtual double foundEvent() const = 0 ;
};

/** */
enum ThresholdEventType {AnyAbove, AllAbove, AnyBelow, AllBelow};

/** */
class ThresholdEventDetector : public EventDetectorBase
{
public:
  /** */
  ThresholdEventDetector(double threshold, ThresholdEventType eventType,
    bool terminateOnDetection=false)
    : threshold_(threshold), eventType_(eventType),
      gotIt_(false), eventTime_(-1.0e300),
      terminateOnDetection_(terminateOnDetection) {}

  /** */
  bool terminateOnDetection() const {return terminateOnDetection_;}

  /** */
  bool checkForEvent(
    const double& t1, const Expr& u1,
    const double& t2, const Expr& u2) ;

  /** */
  double eventTime() const {return eventTime_;}

  /** */
  double foundEvent() const {return gotIt_;}

private:
  double threshold_;
  ThresholdEventType eventType_;
  mutable bool gotIt_;
  mutable double eventTime_;
  bool terminateOnDetection_;
};




}


#endif
