// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "FEApp_Mesh.hpp"

FEApp::Mesh::Mesh() :
  elements()
{
}

FEApp::Mesh::~Mesh()
{
}

void
FEApp::Mesh::addElement(
		   const Teuchos::RefCountPtr<FEApp::AbstractElement>& element)
{
  elements.push_back(element);
}

unsigned int
FEApp::Mesh::numElements() const
{
  return elements.size();
}

FEApp::Mesh::iterator
FEApp::Mesh::begin()
{
  return elements.begin();
}

FEApp::Mesh::const_iterator
FEApp::Mesh::begin() const
{
  return elements.begin();
}

FEApp::Mesh::iterator
FEApp::Mesh::end()
{
  return elements.end();
}

FEApp::Mesh::const_iterator
FEApp::Mesh::end() const
{
  return elements.end();
}
