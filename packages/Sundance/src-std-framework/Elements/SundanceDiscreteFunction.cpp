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

#include "SundanceDiscreteFunction.hpp"
#include "SundanceHomogeneousDOFMap.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceSubtypeEvaluator.hpp"

using namespace SundanceStdMesh;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;


static Time& getLocalValsTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("DF getLocalValues"); 
  return *rtn;
}
static Time& dfCtorTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("DF ctor"); 
  return *rtn;
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
                                   const string& name)
  : DiscreteFunctionStub(name, space.nFunc(),
                         rcp(new DiscreteFunctionData(space))), 
    FuncWithBasis(space.basis()),
    data_()
{
  TimeMonitor timer(dfCtorTimer());
  data_ = rcp_dynamic_cast<DiscreteFunctionData>(dataStub());
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
                                   const Array<string>& name)
  : DiscreteFunctionStub(name, space.nFunc(),
                         rcp(new DiscreteFunctionData(space))), 
    FuncWithBasis(space.basis()),
    data_()
{
  TimeMonitor timer(dfCtorTimer());
  data_ = rcp_dynamic_cast<DiscreteFunctionData>(dataStub());
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
                                   const double& constantValue,
                                   const string& name)
  : DiscreteFunctionStub(name, space.nFunc(),
                         rcp(new DiscreteFunctionData(space, constantValue))), 
    FuncWithBasis(space.basis()),
    data_()
{
  TimeMonitor timer(dfCtorTimer());
  data_ = rcp_dynamic_cast<DiscreteFunctionData>(dataStub());
  Vector<double> vec = data_->getVector();
  vec.setToConstant(constantValue);
  data_->setVector(vec);
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
                                   const double& constantValue,
                                   const Array<string>& name)
  : DiscreteFunctionStub(name, space.nFunc(),
                         rcp(new DiscreteFunctionData(space, constantValue))), 
    FuncWithBasis(space.basis()),
    data_()
{
  TimeMonitor timer(dfCtorTimer());
  data_ = rcp_dynamic_cast<DiscreteFunctionData>(dataStub());
  Vector<double> vec = data_->getVector();
  vec.setToConstant(constantValue);
  data_->setVector(vec);
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
                                   const Vector<double>& vec,
                                   const string& name)
  : DiscreteFunctionStub(name, space.nFunc(),
                         rcp(new DiscreteFunctionData(space, vec))), 
    FuncWithBasis(space.basis()),
    data_()
{
  TimeMonitor timer(dfCtorTimer());
  data_ = rcp_dynamic_cast<DiscreteFunctionData>(dataStub());
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
                                   const Vector<double>& vec,
                                   const Array<string>& name)
  : DiscreteFunctionStub(name, space.nFunc(),
                         rcp(new DiscreteFunctionData(space, vec))), 
    FuncWithBasis(space.basis()),
    data_()
{
  TimeMonitor timer(dfCtorTimer());
  data_ = rcp_dynamic_cast<DiscreteFunctionData>(dataStub());
}

void DiscreteFunction::setVector(const Vector<double>& vec) 
{
  data_->setVector(vec);
}

void DiscreteFunction::updateGhosts() const
{
  data_->updateGhosts();
}


void DiscreteFunction::getLocalValues(int cellDim, 
                                      const Array<int>& cellLID,
                                      Array<Array<double> >& localValues) const 
{
  TimeMonitor timer(getLocalValsTimer());
  data_->getLocalValues(cellDim, cellLID, localValues);
}


const DiscreteFunction* DiscreteFunction::discFunc(const Expr& expr)
{
  const ExprBase* e = expr.ptr().get();
  const DiscreteFunction* df 
    = dynamic_cast<const DiscreteFunction*>(e);

  TEST_FOR_EXCEPTION(df==0, RuntimeError,
                     "failed to cast " << expr << " to a discrete function");

  return df;
}



DiscreteFunction* DiscreteFunction::discFunc(Expr& expr)
{
  DiscreteFunction* df 
    = dynamic_cast<DiscreteFunction*>(expr.ptr().get());

  TEST_FOR_EXCEPTION(df==0, RuntimeError,
                     "failed to cast " << expr << " to a discrete function");

  return df;
}



