/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceCoordDeriv.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

DiscreteFuncElement::DiscreteFuncElement(DiscreteFunctionStub* master, 
                                         const string& name,
                                         const string& suffix,
                                         int myIndex)
	: LeafExpr(), 
    FuncElementBase(name, suffix),
    master_(master),
    myIndex_(myIndex)
{}

void DiscreteFuncElement::findNonzeros(const EvalContext& context,
                                       const Set<MultiIndex>& multiIndices,
                                       const Set<MultiSet<int> >& activeFuncIDs,
                                       const Set<int>& allFuncIDs,
                                       bool regardFuncsAsConstant) const
{

  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for discrete func " 
                       << toString() << " subject to multiindices "
                       << multiIndices);

  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       allFuncIDs, regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }


  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices);
  
  for (Set<MultiIndex>::const_iterator 
         i=multiIndices.begin(); i != multiIndices.end(); i++)
    {
      if (i->order()==1)
        {
          subset->addDeriv(new CoordDeriv(i->firstOrderDirection()), 
                           VectorDeriv);
        }
      if (i->order()==0)
        {
          subset->addDeriv(MultipleDeriv(),
                           VectorDeriv);
        }
    }

  addKnownNonzero(context, multiIndices, activeFuncIDs,
                  allFuncIDs, regardFuncsAsConstant);
}

XMLObject DiscreteFuncElement::toXML() const 
{
	XMLObject rtn("DiscreteFuncElement");
	rtn.addAttribute("name", name());
	return rtn;
}

