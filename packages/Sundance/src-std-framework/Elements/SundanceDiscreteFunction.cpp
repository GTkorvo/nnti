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
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceSubtypeEvaluator.hpp"
#include "PlayaDefaultBlockVectorDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaDefaultBlockVectorImpl.hpp"
#endif

namespace Sundance
{
using namespace Teuchos;
using std::string;
using std::runtime_error;
using std::endl;

static Time& getLocalValsTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("DF getLocalValues"); 
  return *rtn;
}
static Time& dfCtorTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("DF ctor"); 
  return *rtn;
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
  const string& name)
  : DiscreteFunctionStub(tuple(name), space.dimStructure(),
    getRCP(new DiscreteFunctionData(space))), 
    FuncWithBasis(space.basis()),
    data_()
{
  TimeMonitor timer(dfCtorTimer());
  data_ = rcp_dynamic_cast<DiscreteFunctionData>(dataStub());
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
  const Array<string>& name)
  : DiscreteFunctionStub(name, space.dimStructure(),
    getRCP(new DiscreteFunctionData(space))), 
    FuncWithBasis(space.basis()),
    data_()
{
  TimeMonitor timer(dfCtorTimer());
  data_ = rcp_dynamic_cast<DiscreteFunctionData>(dataStub());
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
  const double& constantValue,
  const string& name)
  : DiscreteFunctionStub(tuple(name), space.dimStructure(),
    getRCP(new DiscreteFunctionData(space, constantValue))), 
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
  : DiscreteFunctionStub(name, space.dimStructure(),
    getRCP(new DiscreteFunctionData(space, constantValue))), 
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
  : DiscreteFunctionStub(tuple(name), space.dimStructure(),
    getRCP(new DiscreteFunctionData(space, vec))), 
    FuncWithBasis(space.basis()),
    data_()
{
  TimeMonitor timer(dfCtorTimer());
  data_ = rcp_dynamic_cast<DiscreteFunctionData>(dataStub());
}

DiscreteFunction::DiscreteFunction(const DiscreteSpace& space, 
  const Vector<double>& vec,
  const Array<string>& name)
  : DiscreteFunctionStub(name, space.dimStructure(),
    getRCP(new DiscreteFunctionData(space, vec))), 
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


RCP<const MapStructure> DiscreteFunction::getLocalValues(int cellDim, 
  const Array<int>& cellLID,
  Array<Array<double> >& localValues) const 
{
  TimeMonitor timer(getLocalValsTimer());
  return data_->getLocalValues(cellDim, cellLID, localValues);
}


const DiscreteFunction* DiscreteFunction::discFunc(const Expr& expr)
{
  const ExprBase* e = expr.ptr().get();
  const DiscreteFunction* df 
    = dynamic_cast<const DiscreteFunction*>(e);

  TEUCHOS_TEST_FOR_EXCEPTION(df==0, runtime_error,
    "failed to cast " << expr << " to a discrete function. "
    "It appears to be of type " << e->typeName());

  return df;
}



DiscreteFunction* DiscreteFunction::discFunc(Expr& expr)
{
  DiscreteFunction* df 
    = dynamic_cast<DiscreteFunction*>(expr.ptr().get());

  TEUCHOS_TEST_FOR_EXCEPTION(df==0, runtime_error,
    "failed to cast " << expr << " to a discrete function. "
    "It appears to be of type " << expr.ptr()->typeName());

  return df;
}


RCP<DiscreteFuncDataStub> DiscreteFunction::getRCP(DiscreteFunctionData* ptr)
{
  return rcp_dynamic_cast<DiscreteFuncDataStub>(rcp(ptr));
}





void updateDiscreteFunction(const Expr& newVals, Expr old)
{
  Vector<double> vIn = getDiscreteFunctionVector(newVals);
  Vector<double> vOut = getDiscreteFunctionVector(old);

  vOut.acceptCopyOf(vIn);
  setDiscreteFunctionVector(old, vOut);
}

void addVecToDiscreteFunction(Expr df, const Vector<double>& v)
{
  Vector<double> dfVec = getDiscreteFunctionVector(df);
  dfVec.update(1.0, v);
  setDiscreteFunctionVector(df, dfVec);
}

Expr copyDiscreteFunction(const Expr& u0, const string& name)
{
  const DiscreteFunction* df 
    = dynamic_cast<const DiscreteFunction*>(u0.ptr().get());

  /* Case 1: u is a discrete function */
  if (df != 0)
  {
    Vector<double> dfVec = df->getVector().copy();
    return new DiscreteFunction(df->discreteSpace(), dfVec, name);
  }

  /* Case 2: u is an element of a length-one expression whose DiscreteFunction
   * wrapper has gotten lost in dereferencing. */
  const DiscreteFuncElement* dfe 
    = dynamic_cast<const DiscreteFuncElement*>(u0.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(dfe!=0 && u0.size() > 1, runtime_error,
    "attempt to access vector of a single element of a multicomponent "
    "DiscreteFunction");
  if (dfe != 0)
  {
    Vector<double> dfVec 
      = DiscreteFunctionData::getData(dfe)->getVector().copy();
    return new DiscreteFunction(
      DiscreteFunctionData::getData(dfe)->discreteSpace(),
      dfVec, name);
  }

  /* Case 3: u is a list of discrete functions */
  Array<Expr> rtn(u0.size());
  for (int b=0; b<u0.size(); b++)
  {
    rtn[b] = copyDiscreteFunction(u0[b], 
      name + "[" + Teuchos::toString(b) + "]");
  }
  return new ListExpr(rtn);
}

void setDiscreteFunctionVector(Expr u, const Vector<double>& v)
{
  DiscreteFunction* df 
    = dynamic_cast<DiscreteFunction*>(u.ptr().get());

  /* Case 1: u is a discrete function */
  if (df != 0)
  {
    df->setVector(v);
    return;
  }

  /* Case 2: u is an element of a length-one expression whose DiscreteFunction
   * wrapper has gotten lost in list element dereferencing. */
  DiscreteFuncElement* dfe 
    = dynamic_cast<DiscreteFuncElement*>(u.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(dfe!=0 && u.size() > 1, runtime_error,
    "attempt to set vector of a single element of a multicomponent "
    "DiscreteFunction");
  if (dfe != 0)
  {
    DiscreteFunctionData::getData(dfe)->setVector(v);
    return;
  }

  /* At this point, the vector should be a block vector */
  TEUCHOS_TEST_FOR_EXCEPTION((df==0 && dfe==0) && u.size()==1, 
    runtime_error,
    "non-block vector should be a discrete function in setDFVector()");

  /* Case 3: u is a list of discrete functions */
  for (int b=0; b<u.size(); b++)
  {
    setDiscreteFunctionVector(u[b], v.getBlock(b));
  }
}

Vector<double> getDiscreteFunctionVector(const Expr& u)
{
  const DiscreteFunction* df 
    = dynamic_cast<const DiscreteFunction*>(u.ptr().get());

  /* Case 1: u is a discrete function */
  if (df != 0)
  {
    return df->getVector();
  }

  /* Case 2: u is an element of a length-one expression whose DiscreteFunction
   * wrapper has gotten lost in dereferencing. */
  const DiscreteFuncElement* dfe 
    = dynamic_cast<const DiscreteFuncElement*>(u.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(dfe!=0 && u.size() > 1, runtime_error,
    "attempt to access vector of a single element of a multicomponent "
    "DiscreteFunction");
  if (dfe != 0)
  {
    return DiscreteFunctionData::getData(dfe)->getVector();
  }

  /* Case 3: u is a list of discrete functions */
  TEUCHOS_TEST_FOR_EXCEPTION(df==0 && u.size()==1, runtime_error,
    "non-block vector should be a discrete function in getDiscreteFunctionVector()");
  Array<Vector<double> > vec(u.size());
  for (int b=0; b<u.size(); b++)
  {
    vec[b] = getDiscreteFunctionVector(u[b]);
  }
  return blockVector(vec);
}


Mesh getDiscreteFunctionMesh(const Expr& u)
{
  const DiscreteFunction* df 
    = dynamic_cast<const DiscreteFunction*>(u.ptr().get());

  /* Case 1: u is a discrete function */
  if (df != 0)
  {
    return df->mesh();
  }

  /* Case 2: u is an element of a length-one expression whose DiscreteFunction
   * wrapper has gotten lost in dereferencing. */
  const DiscreteFuncElement* dfe 
    = dynamic_cast<const DiscreteFuncElement*>(u.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(dfe!=0 && u.size() > 1, runtime_error,
    "attempt to access vector of a single element of a multicomponent "
    "DiscreteFunction");
  if (dfe != 0)
  {
    return DiscreteFunctionData::getData(dfe)->mesh();
  }

  /* Case 3: u is a list of discrete functions */
  TEUCHOS_TEST_FOR_EXCEPTION(df==0 && u.size()==1, runtime_error,
    "non-block vector should be a discrete function in getDiscreteFunctionVector()");
  return getDiscreteFunctionMesh(u[0]);
}


}

