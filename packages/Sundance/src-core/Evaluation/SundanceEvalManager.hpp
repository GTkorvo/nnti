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

#ifndef SUNDANCE_EVALMANAGER_H
#define SUNDANCE_EVALMANAGER_H

#include "SundanceDefs.hpp"
#include "SundanceEvalContext.hpp"
#include "SundanceAbstractEvalMediator.hpp"
#include "SundanceTempStack.hpp"
#include "SundanceNoncopyable.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  class CoordExpr;


  namespace Internal
    {
      class MultiIndex;
      class DiscreteFuncElement;

      using namespace Internal;

      /**
       * EvalManager provides methods for interfacing to the framework
       * through an AbstractEvalMediator and managing temporary variables
       * through a TempStack.
       *
       * If no mediator is set, string evaluations will be done 
       */
      class EvalManager : public Noncopyable
        {
        public:
          /** Empty ctor */
          EvalManager();

          /** */
          void evalCoordExpr(const CoordExpr* expr,
                             RefCountPtr<EvalVector>&  result) const ;

          /** */
          void evalCellDiameterExpr(const CellDiameterExpr* expr,
                                    RefCountPtr<EvalVector>&  result) const ;

          /** */
          void evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                                       const Array<MultiIndex>& mi,
                                       Array<RefCountPtr<EvalVector> >& result) const ;

          /** */
          void setMediator(const RefCountPtr<AbstractEvalMediator>& med) 
          {mediator_ = med;}

          /** */
          void setVecSize(int vecSize) {stack().setVecSize(vecSize);}
          

          /** Return a pointer to the mediator. We'll need the
           * mediator for computing framework-specific functions.
           */
          const AbstractEvalMediator* mediator() const {return mediator_.get();}

          /** */
          void setRegion(const EvalContext& region)
            {region_ = region;}

          /** */
          const EvalContext& getRegion() const {return region_;}

          /** */
          static TempStack& stack();

          /** */
          int getMaxDiffOrder() const ;


          /** */
          RefCountPtr<EvalVector> popVector() const ;

        private:

          EvalContext region_;

          RefCountPtr<AbstractEvalMediator> mediator_;

      };

    }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
