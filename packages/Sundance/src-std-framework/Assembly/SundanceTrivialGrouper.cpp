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

#include "SundanceTrivialGrouper.hpp"
#include "SundanceRefIntegral.hpp"
#include "SundanceQuadratureIntegral.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;

void TrivialGrouper::findGroups(const EquationSet& eqn,
                                const CellType& cellType,
                                int cellDim,
                                const QuadratureFamily& quad,
                                const RefCountPtr<SparsitySuperset>& sparsity,
                                Array<IntegralGroup>& groups) const
{
  Tabs tab;
  VerbositySetting verb = GrouperBase::classVerbosity();
  SUNDANCE_OUT(verb > VerbLow, 
               tab << "trivial grouper num derivs = " << sparsity->numDerivs() << endl);
  SUNDANCE_OUT(verb > VerbLow, 
               tab << "cell type = " << cellType);

  SUNDANCE_OUT(verb > VerbMedium,  
               tab << "sparsity = " << endl << *sparsity << endl);

  int vecCount=0;
  int constCount=0;

  for (int i=0; i<sparsity->numDerivs(); i++)
    {
      Tabs tab1;

      const MultipleDeriv& d = sparsity->deriv(i);
      if (d.order()==0) 
        {
          RefCountPtr<ElementIntegral> integral ;
          int resultIndex;
          if (sparsity->isConstant(i))
            {
              integral = rcp(new RefIntegral(cellDim, cellType));
              resultIndex = constCount++;
            }
          else
            {
              integral = rcp(new QuadratureIntegral(cellDim, cellType, quad));
              resultIndex = vecCount++;
            }
          groups.append(IntegralGroup(tuple(integral),
                                      tuple(tuple(resultIndex))));
        }
      else
        {

          BasisFamily testBasis;
          BasisFamily unkBasis;
          MultiIndex miTest;
          MultiIndex miUnk;
          int testID;
          int unkID;
          bool isOneForm;
          Array<int> alpha;
          Array<int> beta;
          extractWeakForm(eqn, d, testBasis, unkBasis, miTest, miUnk, testID, unkID,
                          isOneForm);

          SUNDANCE_OUT(verb > VerbMedium, 
                       tab1 << "test ID: " << testID);

          SUNDANCE_OUT(!isOneForm && verb > VerbMedium,
                       tab1 << "unk funcID: " << unkID << endl);
                   
          SUNDANCE_OUT(verb > VerbMedium, tab1 << "deriv = " << d);
          SUNDANCE_OUT(verb > VerbMedium && sparsity->isConstant(i), 
                       tab1 << "coeff is constant");
          SUNDANCE_OUT(verb > VerbMedium && !sparsity->isConstant(i), 
                       tab1 << "coeff is non-constant");

          RefCountPtr<ElementIntegral> integral;
          int resultIndex;
          if (sparsity->isConstant(i))
            {
              if (isOneForm)
                {
                  if (miTest.order()==1)
                    {
                      alpha = tuple(miTest.firstOrderDirection());
                    }
                  SUNDANCE_OUT(verb > VerbMedium,
                               tab1 << "creating reference integral for one-form");
                  integral = rcp(new RefIntegral(cellDim, cellType,
                                                 testBasis, alpha, 
                                                 miTest.order()));
                }
              else
                {
                  if (miTest.order()==1)
                    {
                      alpha = tuple(miTest.firstOrderDirection());
                    }
                  if (miUnk.order()==1)
                    {
                      beta = tuple(miUnk.firstOrderDirection());
                    }
                  SUNDANCE_OUT(verb > VerbMedium,
                               tab1 << "creating reference integral for two-form");
                  integral = rcp(new RefIntegral(cellDim, cellType,
                                                 testBasis, alpha, miTest.order(),
                                                 unkBasis, beta, miUnk.order()));
                }
              resultIndex = constCount++;
            }
          else
            {
              if (isOneForm)
                {
                  if (miTest.order()==1)
                    {
                      alpha = tuple(miTest.firstOrderDirection());
                    }
                  SUNDANCE_OUT(verb > VerbMedium,
                               tab1 << "creating quadrature integral for two-form");
                  integral = rcp(new QuadratureIntegral(cellDim, cellType,
                                                        testBasis, alpha, 
                                                        miTest.order(), quad));
                }
              else
                {
                  if (miTest.order()==1)
                    {
                      alpha = tuple(miTest.firstOrderDirection());
                    }
                  if (miUnk.order()==1)
                    {
                      beta = tuple(miUnk.firstOrderDirection());
                    }
                  SUNDANCE_OUT(verb > VerbMedium,
                               tab1 << "creating quadrature integral for two-form");
                  integral = rcp(new QuadratureIntegral(cellDim, cellType,
                                                        testBasis, alpha, 
                                                        miTest.order(),
                                                        unkBasis, beta, 
                                                        miUnk.order(), quad));
                }
              resultIndex = vecCount++;
            }

          if (isOneForm)
            {
              groups.append(IntegralGroup(tuple(testID), tuple(integral),
                                          tuple(tuple(resultIndex))));
            }
          else
            {
              groups.append(IntegralGroup(tuple(testID), tuple(unkID),
                                          tuple(integral),
                                          tuple(tuple(resultIndex))));
            }
        }
    }
}

