/* @HEADER@ */
/* @HEADER@ */


#include "SundanceTempStack.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceSparsityPattern.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::FrameworkInterface;
using namespace Teuchos;

TempStack::TempStack(int vecSize)
  : 
  numericalEval_(true),
  vecSize_(vecSize),
  fullVecStack_(),
  trivialVecStack_(),
  vecArrayStack_(),
  trivialVecsAllocated_(0),
  fullVecsAllocated_(0),
  vecArraysAllocated_(0),
  trivialVecsAccessed_(0),
  fullVecsAccessed_(0),
  vecArraysAccessed_(0)
{}

TempStack::TempStack()
  : 
  numericalEval_(false),
  vecSize_(0),
  fullVecStack_(),
  trivialVecStack_(),
  vecArrayStack_(),
  trivialVecsAllocated_(0),
  fullVecsAllocated_(0),
  vecArraysAllocated_(0),
  trivialVecsAccessed_(0),
  fullVecsAccessed_(0),
  vecArraysAccessed_(0)
{}

RefCountPtr<EvalVector> TempStack::popTrivialVector() 
{
  RefCountPtr<EvalVector> rtn;

  if (trivialVecStack_.empty())
    {
      trivialVecsAllocated_++;
      rtn = rcp(new EvalVector());
      if (!numericalEval_) rtn->doStringEval();
    }
  else
    {
      rtn = trivialVecStack_.top();
      trivialVecStack_.pop();
    }
  trivialVecsAccessed_++;
  
  return rtn;
}

RefCountPtr<EvalVector> TempStack::popFullVector() 
{
  RefCountPtr<EvalVector> rtn;

  if (fullVecStack_.empty())
    {
      fullVecsAllocated_++;
      rtn = rcp(new EvalVector(vecSize_));
      if (!numericalEval_) rtn->doStringEval();
    }
  else
    {
      rtn = fullVecStack_.top();
      rtn->resize(vecSize_);
      fullVecStack_.pop();
    }
  fullVecsAccessed_++;
  
  return rtn;
}

void TempStack::pushVector(const RefCountPtr<EvalVector>& vec) 
{
  if (vec->length()==0)
    {
      trivialVecStack_.push(vec);
    }
  else
    {
      fullVecStack_.push(vec);
    }
}

RefCountPtr<EvalVectorArray> 
TempStack::popVectorArray(const SparsityPattern* sparsity) 
{
  RefCountPtr<EvalVectorArray> rtn;

  if (vecArrayStack_.empty())
    {
      vecArraysAllocated_++;
      rtn = rcp(new EvalVectorArray());
    }
  else
    {
      rtn = vecArrayStack_.top();
      vecArrayStack_.pop();
    }

  rtn->resize(sparsity->numDerivs());

  for (int i=0; i<sparsity->numDerivs(); i++)
    {
      if (sparsity->isZero(i) || sparsity->isZero(i))
        {
          (*rtn)[i] = popTrivialVector();
        }
      else
        {
          (*rtn)[i] = popFullVector();
        }
    }
  
  vecArraysAccessed_++;

  return rtn;
}

void TempStack::pushVectorArray(const RefCountPtr<EvalVectorArray>& vecs)  
{
  for (int i=0; i<vecs->size(); i++)
    {
      pushVector((*vecs)[i]);
    }
  vecArrayStack_.push(vecs);
}

void TempStack::resetCounters()
{
  trivialVecsAllocated_=0;
  fullVecsAllocated_=0;
  vecArraysAllocated_=0;
  trivialVecsAccessed_=0;
  fullVecsAccessed_=0;
  vecArraysAccessed_=0;
}


