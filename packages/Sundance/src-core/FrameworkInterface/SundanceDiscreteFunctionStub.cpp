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

#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceSpectralExpr.hpp"


using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

DiscreteFunctionStub::DiscreteFunctionStub(const string& name, int nElems,
                                           const RefCountPtr<DiscreteFuncDataStub>& data)
	: ListExpr(), data_(data)
{
  for (int i=0; i<nElems; i++)
    {
      string suffix;
      if (nElems > 1) suffix = "[" + Teuchos::toString(i) + "]";
      append(new DiscreteFuncElement(data, name, suffix, i));
    }
}



DiscreteFunctionStub::DiscreteFunctionStub(const string& name, const SpectralBasis& sbasis, int nElems,
                                           const RefCountPtr<DiscreteFuncDataStub>& data)
  : ListExpr(), data_(data)
{
  int counter = 0;
  for (int i=0; i<nElems; i++)
    {
      string suffix;
      Array<Expr> Coeffs(sbasis.nterms());

      for(int n=0; n< sbasis.nterms(); n++)
        {
          suffix = "[" + Teuchos::toString(counter) + "]";
          Coeffs[n] = new DiscreteFuncElement(data, name, suffix, counter);
          counter++;
        }

      append(new SpectralExpr(sbasis, Coeffs));
    }
}


DiscreteFunctionStub::DiscreteFunctionStub(const Array<string>& names, int nElems,
                                           const RefCountPtr<DiscreteFuncDataStub>& data)
	: ListExpr(), data_(data)
{
  TEST_FOR_EXCEPTION(names.size() != nElems,
                     RuntimeError,
                     "mismatch between size of names array=" << names 
                     << " and specified size=" << nElems);

  for (int i=0; i<nElems; i++)
    {
      append(new DiscreteFuncElement(data, "", names[i], i));
    }
}


DiscreteFunctionStub::DiscreteFunctionStub(const Array<string>& names, 
                                           const SpectralBasis& sbasis, int nElems,
                                           const RefCountPtr<DiscreteFuncDataStub>& data)
  : ListExpr(), data_(data)
{

  int elemCounter = 0;
  for (int i=0; i<nElems; i++)
    {
      TEST_FOR_EXCEPTION(names.size() != nElems,
                         RuntimeError,
                         "mismatch between size of names array=" << names 
                         << " and specified size=" << nElems);
      
      string suffix;
      Array<Expr> Coeffs(sbasis.nterms());
      
      for(int n=0; n< sbasis.nterms(); n++)
        {
          suffix = "[" + Teuchos::toString(n) + "]";
          Coeffs[n] = new DiscreteFuncElement(data, names[i], suffix, elemCounter);
          elemCounter++;
        }
            
      append(new SpectralExpr(sbasis, Coeffs));
    }
}
