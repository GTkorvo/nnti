/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/

#include "TSFEpetraVectorSpace.hpp"
#include "TSFEpetraVector.hpp"
#include "Teuchos_Utils.hpp"

using namespace TSFExtended;
using namespace Teuchos;


EpetraVectorSpace::EpetraVectorSpace()
	: TSFCore::EpetraVectorSpace()
{}

EpetraVectorSpace::EpetraVectorSpace(const RefCountPtr<const Epetra_Map>& localMap)
	: TSFCore::EpetraVectorSpace(localMap)
{}

RefCountPtr<TSFCore::Vector<double> > EpetraVectorSpace::createMember() const
// Vector<double>  EpetraVectorSpace::createMember() const
{
  Epetra_Vector* v = new Epetra_Vector(*epetra_map(), false);

  RefCountPtr<Epetra_Vector> vec = rcp(v);
  //  RefCountPtr<const TSFCore::EpetraVectorSpace> me = rcp(this, false);
  RefCountPtr<const EpetraVectorSpace> me = rcp(this, false);
  return rcp(new EpetraVector(vec, me));
}

// string EpetraVectorSpace::describe() const 
// {
//   return describe(0);
// }

// 	string rtn = "EpetraVectorSpace[";
//   rtn += "nLocal=" 
//     + Teuchos::toString(epetra_map()->NumMyElements())
//     + " nGlobal=" 
//     + Teuchos::toString(epetra_map()->NumGlobalElements()) 
//     + "]";


//   return rtn;
// }


