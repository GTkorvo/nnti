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

#ifndef SUNDANCE_UNARYMINUS_H
#define SUNDANCE_UNARYMINUS_H

#include "SundanceDefs.hpp"
#include "SundanceUnaryExpr.hpp"
#include "SundanceUnaryMinusEvaluator.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using std::string;

  namespace Internal
    {
      /** */
      class UnaryMinus : public UnaryExpr,
                         GenericEvaluatorFactory<UnaryMinus, UnaryMinusEvaluator>
        {
        public:
          /** construct with the argument */
          UnaryMinus(const RefCountPtr<ScalarExpr>& arg);

          /** virtual dtor */
          virtual ~UnaryMinus() {;}


          /** */
          virtual ostream& toText(ostream& os, bool paren) const ;

          /** */
          virtual ostream& toLatex(ostream& os, bool paren) const ;

          /** */
          virtual XMLObject toXML() const ;

          /** */
          virtual bool isLinear() const {return true;}

          /** */
          virtual Set<MultiSet<int> > internalFindQ_W(int order, 
                                                   const EvalContext& context) const ;

          /** */
          virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}

        };
    }
}
                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  



#endif
